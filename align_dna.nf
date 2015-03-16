#!/usr/bin/env nextflow

import CopyHelper
import ChannelUtil
import FastQC
import PathUtil
import ParamsHelper

// ---------------------------------------------------------------------------
// Read preprocessing and alignment for DNA (WES or WGS) reads.
// ---------------------------------------------------------------------------

if (params.verbose)
    echo true

ParamsHelper.checkNonEmptyParam(params.runID, "runID");
ParamsHelper.checkNonEmptyParam(params.runPlatform, "runPlatform");
ParamsHelper.checkNonEmptyParam(params.dataDir, "dataDir");

copyHelper = new CopyHelper(params.dataDir, params.printCopyMsgs)

// Open channel for left and right files and merge it into triples, the
// first entry is the LCS of the file names that can be used as a read
// pair identifier.
readPairs = ChannelUtil.createFilePairChannel(
        params.runID,
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_{R,}1.fastq.gz'].join(File.separator)),
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_{R,}2.fastq.gz'].join(File.separator)),
        )

// Genome and index files.
indexFileBWA = file(params.indexBWA)

// Duplicate the read pairs into one queue for runFastQCOriginal
// and runTrimming.
(readPairsFastQCOriginal, readPairsRunTrimming) = readPairs.separate(2) { x -> [x, x] }

// --------------------------------------------------------------------------
// Step 1a) Run FastQC
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCOriginal {
    cpus params.fastqc.cpus
    module 'fastqc/0.11.2'

    input:
    set runID, file(readL), file(readR) from readPairsFastQCOriginal

    output:
    set file('*.zip'), file('*.html') into fastqcOutputOriginal

    script:
    """
    set -x
    fastqc -t params.fastqc.cpus -o . ${readL} ${readR}
    """
}

copyHelper.copyFiles(fastqcOutputOriginal, 'reports/fastqc-original');

System.exit(0);

// --------------------------------------------------------------------------
// Step 1b) Run adapter trimming
//
// - yields trimmed read, used as downstream input
// --------------------------------------------------------------------------

process runTrimming {
    cpus params.skewer.cpus
    module 'skewer/0.1.124'

    input:
    set runID, file(readL), file(readR) from readPairsRunTrimming

    output:
    set runID, file { "out/${readL}" }, file { "out/${readR}" } into readPairsTrimmed
    set file("*.log") into trimmingLogs

    script:
    """
    set -x
    # call Skewer
    skewer \\
      -x ${params.skewer.adaptersR1} \\
      -y ${params.skewer.adaptersR2} \\
      -m pe \\
      -z \\
      -t ${params.skewer.cpus} \\
      ${readL} \\
      ${readR}
    # compute name of left/right Skewer result file
    NAMEBASE=${readL}
    LEFT=\${NAMEBASE%.gz}-trimmed-pair1.fastq.gz
    RIGHT=\${NAMEBASE%.gz}-trimmed-pair2.fastq.gz
    # move Skewer output to expected file names
    mkdir -p out
    mv \${LEFT} out/${readL}
    mv \${RIGHT} out/${readR}
    """
}

// Duplicate the read pairs into multiple queues for processing / copying out.
(readPairsFastQCTrimmed,
 readPairsRunMapping,
 readPairsTrimmedCopyOut) = readPairsTrimmed.separate(3) { x -> [ x, x, x ] }

// Copy out results from trimming step (map removes the pair).
copyHelper.copyFiles(trimmingLogs, 'reports/trimming');
copyHelper.copyFiles(readPairsTrimmedCopyOut.map { [it[1], it[2]] }, 'fastq/trimmed');

// --------------------------------------------------------------------------
// Step 2a) Run FastQC on trimmed
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCTrimmed {
    cpus params.fastqc.cpus
    module 'fastqc/0.11.2'

    input:
    set runID, file(readL), file(readR) from readPairsFastQCTrimmed

    output:
    set file('*.zip'), file('*.html') into fastqcOutputTrimmed

    script:
    """
    set -x
    fastqc -t ${params.fastqc.cpus} -o . ${readL} ${readR}
    """
}

copyHelper.copyFiles(fastqcOutputTrimmed, 'reports/fastqc-trimmed');

// --------------------------------------------------------------------------
// Step 2b) Align reads using BWA-MEM
//
// - align reads
// - sort
// - mark duplicates
// - yields alignment for downstream processing
// --------------------------------------------------------------------------

// Group trimmed read FASTQ files by runID (the runID is part of the output
// of a previous process).
jointBams = readPairsRunMapping.map{f -> [f[0], f[1], f[2]] }.groupTuple()

// The alignments are written to the temporary files alignment.bam. These
// BAM files are already sorted.
process runReadMapping {
    cpus params.bwa.cpus
    module 'bwa/0.7.12'
    module 'samtools/1.2'
    module 'samblaster/0.1.21'

    input:
    indexFileBWA
    set runID, readL, readR from jointBams

    output:
    file { "${runID}.bam*" } into bamFilesOut
    set runID, file { "${runID}.bam" }, file { "${runID}.bam.bai" } into bamFiles

    script:
    """
    set -x
    bwa mem \\
        -R '@RG\tID:${runID}\tSM:${runID}\tPL:${params.runPlatform}' \\
        -t ${params.bwa.cpus} \\
        ${indexFileBWA} \\
        <(zcat ${readL.join(" ")}) \\
        <(zcat ${readR.join(" ")}) \\
        | samblaster \\
        | samtools view -u -Sb - \\
        | samtools sort - ${runID}
    samtools index ${runID}.bam
    """
}

(bamFilesForCoverage,
 bamFilesForQualimap,
 bamFilesForVariantCalling) = bamFiles.separate(3) { x -> [ x, x, x ] }

copyHelper.copyFiles(bamFilesOut, 'bam')
