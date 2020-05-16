def helpMessage() {
    log.info"""
    # Nanopore genome polishing
    A pipeline for polishing genomes assembled from Oxford Nanopore reads using Racon and Medaka.

    ## Examples
    nextflow run nanopore_polishing.nf \
    --genomes "04-canu-assembly/*.fasta" \
    --trimmedReads "03-trimmed-fastq/*.fastq.gz"
    

    ## Parameters
    --genomes <glob>
        Required
        A glob of the fasta genomes to be polished.
        The basename of the file is used as the genome name.

    --trimmedReads <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --outdir <path>
        Default: `results_nanopore_polishing`
        The directory to store the results in.--g

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.trimmedReads = false
params.genomes = false
params.outdir = "results_nanopore_polishing"

if ( params.genomes ) {
    genomes = Channel
    .fromPath(params.genomes, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { genomesForPolishing }
} else {
    log.info "No genomes supplied, did you include '*.fasta'?"
    exit 1
}

if ( params.trimmedReads ) {
    trimmedReads = Channel
    .fromPath(params.trimmedReads, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { trimmedReadsForPolishing }
} else {
    log.info "No trimmed reads supplied, did you include '*.fastq.gz'?"
    exit 1
}

genomesForPolishing
.combine(trimmedReadsForPolishing, by: 0)
.set { racon }

process versions {
    publishDir "${params.outdir}/", mode: 'copy'

    output:
    file 'versions.txt'

    """
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    echo racon: >> versions.txt
    racon --version >> versions.txt
    echo --------------- >> versions.txt
    echo medaka: >> versions.txt
    medaka --version >> versions.txt
    """

}

// racon parameters as suggested by medaka authors https://github.com/nanoporetech/medaka
process racon {

    publishDir "${params.outdir}/05-racon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from racon

    output:
    set sampleID, "${sampleID}.contigs.racon.fasta", 'input.fastq.gz' into medaka

    """
    minimap2 \
    input.fasta \
    input.fastq.gz > minimap.racon.paf

    racon -m 8 -x -6 -g -8 -w 500 \
    input.fastq.gz \
    minimap.racon.paf \
    input.fasta > ${sampleID}.contigs.racon.fasta

    """
}

// polishing step 2
process medaka {

    publishDir "${params.outdir}/06-medaka-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from medaka

    output:
    set sampleID, "${sampleID}.contigs.racon.medaka.fasta", 'input.fastq.gz' into unknown

    """
    medaka_consensus \
    -d input.fasta \
    -i input.fastq.gz \
    -o ${sampleID}_medaka_output \
    -m r941_min_high_g360

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}.contigs.racon.medaka.fasta

    """
}

log.info "If nothing happened, did you inlcude '*.fasta' and '*.fastq.gz' in the --genomes and --trimmedReads options?"