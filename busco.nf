def helpMessage() {
    log.info"""
    # Busco completeness pipeline
    A pipeline to determine BUSCO completeness using BUSCO V4

    ## Examples
    nextflow run nanopore_polishing.nf \
    --genomes "genomes/*.fasta" \
    --database "ascomycota_odb10"
    

    ## Parameters
    --genomes <glob>
        Required
        A glob of the fasta genomes to be polished.
        The basename of the file is used as the genome name.

    --database <string>
        Optional
        Default is "ascomycota_odb10"

    --outdir <path>
        Default: `results_busco`
        The directory to store the results in.

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false
params.outdir = "results_busco"
params.database = "ascomycota_odb10"

if ( params.genomes ) {
    genomes = Channel
    .fromPath(params.genomes, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { genomesForBusco }
} else {
    log.info "No genomes supplied."
    exit 1
}

// Busco doesn't seem to run if the input file is on rdrive, therefore I first create a local copy 'test.fasta'

process busco {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    set id, "genome.fasta" from genomesForBusco

    output:
    """
    cp genome.fasta ./test.fasta
    docker run -v \$(pwd):/busco_wd ezlabgva/busco:v4.0.5_cv1 busco \
    -i test.fasta \
    -o ${id} \
    -l ${params.database} \
    -m genome \
    -c 16 
    cp ${id}/short_summary*.txt .
    """
}

busco_output
.collect()
.set {all_busco_outputs}

process plotBuscoSummaries {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file 'short_summary.txt' from all_busco_outputs

    output:
    file 'busco_figure.png'
    
    """
    mkdir input_derp
    cp ${params.outdir}/*.txt input_derp/
    docker run -v \$(pwd)/input_derp:/busco_wd ezlabgva/busco:v4.0.5_cv1 generate_plot.py -wd .
    cp input_derp/busco_figure.png .
    """
}
