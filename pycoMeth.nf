def helpMessage() {
    log.info"""
    # Nanopore genome polishing
    A pipeline for polishing genomes assembled from Oxford Nanopore reads using Racon and Medaka.

    ## Examples
    nextflow run nanopore_polishing.nf \
    --genomes "04-canu-assembly/*.fasta" \
    --trimmedReads "03-trimmed-fastq/*.fastq.gz"
    

    ## Parameters
    --reference <glob>
        Required
        A glob of the fasta genomes to be polished.
        The basename of the file is used as the genome name.

    --nanopolish <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --outdir <path>
        Default: `results_nanopore_polishing`
        The directory to store the results in.

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


params.nanopolish = false
params.reference = false
params.outdir = "results_nanopore_polishing"

if ( params.reference ) {
    reference = Channel
    .fromPath(params.reference, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{referenceGenome}
} else {
    log.info "No reference genome supplied, did you include '*.fasta'?"
    exit 1
}

if ( params.nanopolish ) {
    nanopolish = Channel
    .fromPath(params.nanopolish, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{nanopolishOutput}
} else {
    log.info "No nanopolish reads supplied, did you include '*.tsv?"
    exit 1
}

referenceGenome.println
nanopolishOutput.println
