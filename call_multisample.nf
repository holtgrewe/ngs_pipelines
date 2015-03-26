#!/usr/bin/env nextflow

import CopyHelper
import ChannelUtil
import FastQC
import PathUtil
import ParamsHelper

// ---------------------------------------------------------------------------
// Variant calling from multiple samples.
// ---------------------------------------------------------------------------

// Notes: variant-calling of RNA-seq data does not work with Freebayes

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

// Get handle to BED files.
bedExonsUCSCFile = file(params.exonsUCSC)
bedExonsCCDSFile = file(params.exonsCCDS)

// ---------------------------------------------------------------------------
// Perform variant calling with Freebayes
// ---------------------------------------------------------------------------

process runFreebayes {
    cpus params.freebayes.cpus
    module 'freebayes/0.9.20'
    module 'htslib/1.2.1'
    module 'vcflib/b1e9b3'

    input:
    genomeFile
    set poolID, bamFiles from bamFilesIn

    output:
    file { "${params.poolID}.vcf.gz"} into vcfFilesForCutting
    file { "${params.poolID}.vcf.gz*" } into vcfFilesOut

    script:
    """
    #!/bin/bash
    set -x

    # Run freebayes on a certain region on the genome.
    #
    # The --min-base-quality and --min-mapping-quality settings come from the
    # Freebayes author https://www.biostars.org/p/85400/#87189
    run_freebayes()
    {
        set -x
        REGION=\$1
        shift
        freebayes \\
            --genotype-qualities \\
            --min-base-quality 3 \\
            --min-mapping-quality 1 \\
            -r \${REGION} \\
            -f ${genomeFile} \\
            ${bamFiles.join(" ")}
    }
    export -f run_freebayes

    # The length of the windows to use for parallelization.
    WINDOW_LEN=100000

    /usr/bin/time \\
        parallel -t -k -j ${params.freebayes.cpus} run_freebayes :::: \\
        <(fasta_generate_regions.py ${genomeFile}.fai \${WINDOW_LEN} || true) \\
        | vcffirstheader \\
        | bgzip -c /dev/stdin \\
        > ${poolID}.vcf.gz

    # Index the resulting file using tabix
    tabix ${poolID}.vcf.gz
    """
}

copyHelper.copyFiles(vcfFilesOut, 'vcf');

// ---------------------------------------------------------------------------
// Cut down the VCF files to exons
// ---------------------------------------------------------------------------

process cutVCFToExons {
    module 'bedtools2/2.23.0'

    input:
    file vcfFile from vcfFilesForCutting
    bedExonsUCSCFile
    bedExonsCCDSFile

    output:
    file { "out.d/*" } into cutVCFFiles

    script:
    """
    set -x
    mkdir -p out.d

    FILENAME=\$(basename ${vcfFile} .vcf.gz)

    bedtools \\
        intersect \\
        -header \\
        -a ${vcfFile} \\
        -b ${bedExonsUCSCFile} \\
        | bgzip -c /dev/stdin \\
        > out.d/\${FILENAME}.UCSC.vcf.gz

    bedtools \\
        intersect \\
        -header \\
        -a ${vcfFile} \\
        -b ${bedExonsCCDSFile} \\
        | bgzip -c /dev/stdin \\
        > out.d/\${FILENAME}.CCDS.vcf.gz

    tabix out.d/\${FILENAME}.UCSC.vcf.gz
    tabix out.d/\${FILENAME}.CCDS.vcf.gz
    """
}

copyHelper.copyFiles(cutVCFFiles, 'vcf');
