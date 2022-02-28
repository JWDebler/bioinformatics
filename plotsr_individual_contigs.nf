def helpMessage() {
    log.info"""
    # Gene annotaion with GeneMark-Es
    A pipeline to do de novo gene annotation with GeneMark-ES

    ## Examples
    nextflow run annotate_genomes.nf \
    --genomes "genomes/*.fasta"
    

    ## Parameters
    --genomes <glob>
        Required
        A glob of the fasta genomes to be annotated.
        The basename of the file is used as the genome name.

    --cores <int>
        Optional
        Default: 14

    --outdir <path>
        Default: `results_genemark`
        The directory to store the results in.

    --species <string>
        Default "Alentis"
        Species to use for Augustus

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genome1 = false
params.genome2 = false
params.outdir = "plotsr_by_contig"
params.cores = 14

if ( params.genome1 ) {
    genomes = Channel
    .fromPath(params.genomes, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{reference}
} else {
    log.info "No reference supplied."
    exit 1
}

if ( params.genome2 ) {
    genomes = Channel
    .fromPath(params.genomes, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{query}
} else {
    log.info "No query supplied."
    exit 1
}

reference
.merge(query)
.set {reference_and_query}
.view()

return

process extract_contigs {
    tag {sampleID}

    input:
    set sampleID, 'genome.fasta' from reference_and_query

    output:
    set sampleID, "${sampleID}.fasta" into contigsForGenemark, contigsForAugustus //dos2unixed

    """
    dos2unix -n genome.fasta ${sampleID}.fasta
    """
}
