def helpMessage() {
    log.info"""
    # Busco completeness pipeline
    A pipeline to determine BUSCO completeness using BUSCO V4.
    It uses the simple filename (everything in front of the first '.') of your input data for plots, so make sure it is unique.

    ## Examples
    nextflow run busco.nf \
    --genomes "genomes/*.fasta" \
    --database "ascomycota_odb10"
    

    ## Parameters
    --genomes <glob>
        Required
        A glob of the fasta genomes to be polished.
        The simple filename of the file is used as the genome name.

    --database <string>
        Optional
        Default is "ascomycota_odb10"

    --mode <string>
        Optional
        Default: genome
        Available: transcriptome, proteins

    --cores <int>
        Optional
        Default: 16

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
params.mode = "genome"
params.cores = 16

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
    set sampleID, "${sampleID}.fasta" from genomesForBusco

    output:
    file "short_summary*.txt" into busco_output
    """
    cp ${sampleID}.fasta ./${sampleID}.local.fasta
    docker run -v \$(pwd):/busco_wd ezlabgva/busco:v5.2.2_cv1 busco \
    -i ${sampleID}.local.fasta \
    -o ${sampleID} \
    -l ${params.database} \
    -m ${params.mode} \
    -c ${params.cores} 
    cp ${sampleID}/short_summary*.txt .
    """
}

busco_output
.collect()
.set {all_busco_outputs}

process plotBuscoSummaries {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file all_busco_outputs

    output:
    file 'busco_figure.png'
    
    """
    mkdir input_derp
    cp $all_busco_outputs input_derp/
    docker run -v \$(pwd)/input_derp:/busco_wd ezlabgva/busco:v4.0.5_cv1 generate_plot.py -wd .
    cp input_derp/busco_figure.png .
    """
}l
