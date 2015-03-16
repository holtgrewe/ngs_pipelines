#!/usr/bin/env nextflow

import CopyHelper
import ChannelUtil
import FastQC
import PathUtil
import ParamsHelper

// ---------------------------------------------------------------------------
// Variant calling from multiple samples.
// ---------------------------------------------------------------------------

if (params.verbose)
    echo true

ParamsHelper.checkNonEmptyParam(params.inputBam, "inputBam");  // separate by colon ':'
ParamsHelper.checkNonEmptyParam(params.dataDir, "dataDir");
ParamsHelper.checkNonEmptyParam(params.dataDir, "poolID");

copyHelper = new CopyHelper(params.dataDir, params.printCopyMsgs)

// Open channel from input BAM paths and convert to files.
bamFilesIn = Channel.from(params.inputBam.split(':')).map { [params.poolID, file(it)] }.groupTuple()

// Get handle to genome file.
genomeFile = file(params.genome)

// ---------------------------------------------------------------------------
// Perform variant calling with Freebayes
// ---------------------------------------------------------------------------

process runFastQCOriginal {
    cpus params.fastqc.cpus
    module 'freebayes/0.9.21'
    module 'htslib/1.2.1'

    input:
    genomeFile
    set poolID, bamFiles from bamFilesIn

    output:
    file { "${params.poolID}.vcf*" } into vcfFilesOut

    script:
    """
    set -x
    # The --min-base-quality and --min-mapping-quality settings come from the
    # Freebayes author https://www.biostars.org/p/85400/#87189
    freebayes \\
        --genotype-qualities \\
        --min-base-quality 3 \\
        --min-mapping-quality 1 \\
        -f ${genomeFile} \\
        ${bamFiles.join(" ")} \\
        | bgzip -c /dev/stdin \\
        > ${poolID}.vcf.gz
    tabix ${poolID}.vcf.gz
    """
}
