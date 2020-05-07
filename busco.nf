params.workdir = "/home/ubuntu/2020-04-02_nanopore_basecalling/06-medaka-polish"
params.input = "${params.workdir}/*.fasta"
params.outdir = "${params.workdir}/BUSCO_summaries"

sequences = Channel
.fromPath(params.input)
.map{[it.getSimpleName(), it]}

// Busco doesn't seem to run if the input file is on rdrive, therefore I first create a local copy 'test.fasta'

process busco {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    set id, "genome.fasta" from sequences

    output:
    """
    cp genome.fasta ./test.fasta
    docker run -v \$(pwd):/busco_wd ezlabgva/busco:v4.0.5_cv1 busco \
    -i test.fasta \
    -o ${id} \
    -l ascomycota_odb10 \
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
