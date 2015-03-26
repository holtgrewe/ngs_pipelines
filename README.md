NGS Pipelines
=============

Getting Started
---------------

First, install the dependencies:

* Nextflow (http://nextflow.io)
* FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Skewer (https://github.com/relipmoc/skewer)
* BWA (https://github.com/lh3/bwa)
* Samblaster (https://github.com/lh3/bwa)
* Samtools (http://www.htslib.org/)
* STAR aligner (https://github.com/alexdobin/STAR)
* Freebayes (https://github.com/ekg/freebayes)
* htslib (http://www.htslib.org/)
* vcflib (https://github.com/ekg/vcflib)
* bedtools2 (https://github.com/arq5x/bedtools2)
* GNU parallel (should have a package in your Linux distribution)

Then, clone the repository:

    # git clone https://github.com/bihealth/ngs_pipelines.git

And run the pipelines

    # nextflow run align_dna.nf \
        --dataDir $HOME/Data/2015_03_11_triple_cohort/wes_tumor \
        --runID MM065_wes_tumor \
        --runPlatform Illumina \
        -resume \
        -with-trace
    # nextflow run align_dna.nf \
        --dataDir $HOME/Data/2015_03_11_triple_cohort/wes_blood \
        --runID MM065_wes_blood \
        --runPlatform Illumina \
        -resume \
        -with-trace
    # nextflow run align_rna.nf \
        --runID MM065_rna_tumor \
        --runPlatform Illumina \
        --dataDir $HOME/Data/2015_03_11_triple_cohort/rna_tumor \
        -resume \
        -with-trace
    # nextflow run call_multisample.nf \
        --dataDir $HOME/Data/2015_03_11_triple_cohort/variant_calling \
        --inputBam $HOME/Data/2015_03_11_triple_cohort/wes_tumor/bam/MM065_wes_tumor.bam:$HOME/Data/2015_03_11_triple_cohort/wes_blood/bam/MM065_wes_blood.bam \
        --poolID MM065_wes_tumor.MM065_wes_blood \
        -resume \
        -with-trace

Directory Structure
-------------------

The pipeline expects to have one project directory (`2015_03_11_triple_cohort` in the example above) with one sub-directory for each sequenced sample (`wes_blood`, `wes_tumor`, and `rna_tumor`).

For the **input**, each sample folder should have a subdirectory `fastq/original` in which the original FASTQ files reside. Currently, only paired reads are supported, the left reads ("first read in pair") should have a name matching the pattern `*_1.fastq.gz` or `*_R1.fastq.gz`. The second read should then have the name `${NAME}_2.fastq.gz` or `${NAME}_R2.fastq.gz` where `${NAME}` is the prefix of the first read. There can be multiple read pairs in the input directory.

    SAMPLE_OR_WETLAB_ID
    `-- fastq
        `-- original

After running both the alignment and the variant calling, the resulting folder structure will look as follows:

    2015_03_11_triple_cohort
    |-- rna_tumor
    |   |-- bam
    |   |-- fastq
    |   |   |-- original
    |   |   `-- trimmed
    |   `-- reports
    |       |-- alignment
    |       |-- fastqc-original
    |       |-- fastqc-trimmed
    |       |-- split_at_n
    |       `-- trimming
    |-- variant_calling
    |   `-- vcf
    |-- wes_blood
    |   |-- bam
    |   |-- fastq
    |   |   |-- original
    |   |   `-- trimmed
    |   `-- reports
    |       |-- fastqc-original
    |       |-- fastqc-trimmed
    |       `-- trimming
    `-- wes_tumor
        |-- bam
        |-- fastq
        |   |-- original
        |   `-- trimmed
        `-- reports
            |-- fastqc-original
            |-- fastqc-trimmed
            `-- trimming

For each sample, the following subdirectories exist:

* `fastq/trimmed` FASTQ files after adapter trimming.
* `bam` Aligned reads in BAM format.
* `reports` Logs from the trimming/alignment and QC reports.

Also, there is the folder `variant_calling` in the project main directory that contains a `vcf` folder that has the overall calls and calls filtered down to the UCSC and CCDS exons.
