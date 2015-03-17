#!/usr/bin/env nextflow

import CopyHelper
import ChannelUtil
import FastQC
import PathUtil
import ParamsHelper

// ---------------------------------------------------------------------------
// Read preprocessing and alignment for RNA-seq reads
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
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_1.fastq.gz'].join(File.separator)),
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_2.fastq.gz'].join(File.separator)),
        )

// Genome and index files.
genomeFile = file(params.genome)
indexFileSTAR = file(params.indexSTAR)

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
    cpus params.star.cpus
    module 'star/2.4.0j'
    module 'samtools/1.2'
    module 'samblaster/0.1.21'

    input:
    indexFileSTAR
    set runID, readL, readR from jointBams

    output:
    file { "${runID}.bam*" } into bamFilesOut
    set file('out.d/Log.*') into starLogFiles
    set runID, file { "${runID}.bam" } into bamFilesPreSplitting

    script:
    """
    set -x
    mkdir -p out.d
    # Run the STAR RNA-seq aligner.
    #
    # Using the default ENCODE standard settings from the manual below
    # (--outFilterType..--alignMatesGapMax). Also, using Cufflink/Cuffdiff
    # default options (--outSAMstrandField..--outFilterIntronMotifs).
    # Further, we add read groups (--outSAMattrRGline).
    STAR \\
        --genomeDir ${indexFileSTAR} \\
        --runThreadN ${params.star.cpus} \\
        --readFilesIn \\
            <(zcat ${readL.join(" ")}) \\
            <(zcat ${readR.join(" ")}) \\
        --outFileNamePrefix out.d/ \\
        --outFilterType BySJout \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --outSAMattrRGline "ID:${runID}" "SM:${runID}" "PL:${params.runPlatform}"
    # Mask duplicates, sort, convert to BAM.
    samblaster \\
        -i out.d/Aligned.out.sam \\
        | samtools view -Sbu - \\
        | samtools sort -@ ${params.star.cpus} - ${runID}
    # Index resulting BAM file
    samtools index ${runID}.bam
    """
}

copyHelper.copyFiles(bamFilesOut, 'bam')
copyHelper.copyFiles(starLogFiles, 'reports/alignment')

/** DOES NOT HELP CALLING WITH FREEBAYES
// --------------------------------------------------------------------------
// Step 3) Split reads at N in CIGAR string
// --------------------------------------------------------------------------

// For RNA-seq, it is worth performing a splitting at the Ns in CIGAR reads
// after the alignment as described in [1].  Base-callers such as Freebayes
// get confused otherwise.
//
// [1] https://www.broadinstitute.org/gatk/guide/article?id=3891

process runReadSplittingAtNs {
    module 'gatk/3.3-0'

    input:
    genomeFile
    set runID, bamFile from bamFilesPreSplitting

    output:
    file { "${runID}.splitAtN.bam*" } into bamFilesAfterSplitting
    file { 'log/*' } into splitAtNLogFiles

    script:
    """
    mkdir log
    java org.broadinstitute.gatk.engine.CommandLineGATK \\
        -T SplitNCigarReads \\
        -R ${genomeFile} \\
        -I ${bamFile} \\
        -o ${runID}.splitAtN.bam \\
        -rf ReassignOneMappingQuality \\
        -RMQF 255 \\
        -RMQT 60 \\
        -U ALLOW_N_CIGAR_READS \\
        > log/SplitNCigarReads.stdout \\
        2> log/SplitNCigarReads.stderr
    mv ${runID}.splitAtN.bai ${runID}.splitAtN.bam.bai
    """
}

copyHelper.copyFiles(bamFilesAfterSplitting, 'bam')
copyHelper.copyFiles(splitAtNLogFiles, 'reports/split_at_n')
*/
