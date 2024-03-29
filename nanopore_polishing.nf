// run from outside conda environment! If medaka fails it medaka

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
        The basename of the genome and the reads needs to be the same.
        e.g. 'genome_x.fasta' and 'genome_x.fastq.gz'

    --trimmedReads <glob>
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
.set { dos2unix }

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
    """

}

// to handle line end characteres created by exporting fasta from geneious
process dos2unix {
    input:
    set sampleID, "${sampleID}.fasta", "${sampleID}.fastq.gz" from dos2unix

    output:
    set sampleID, "${sampleID}.dos2unix.fasta", "${sampleID}.fastq.gz" into racon

    """
    dos2unix -n ${sampleID}.fasta ${sampleID}.dos2unix.fasta
    """
}


// racon parameters as suggested by medaka authors https://github.com/nanoporetech/medaka
process racon {

    publishDir "${params.outdir}/05-racon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, "${sampleID}.fasta", "${sampleID}.fastq.gz" from racon

    output:
    set sampleID, "${sampleID}.contigs.racon.fasta", "${sampleID}.fastq.gz" into medaka

    """
    minimap2 \
    ${sampleID}.fasta \
    ${sampleID}.fastq.gz > minimap.racon.paf

    racon -m 8 -x -6 -g -8 -w 500 -t 14\
    ${sampleID}.fastq.gz \
    minimap.racon.paf \
    ${sampleID}.fasta > ${sampleID}.contigs.racon.fasta

    """
}

// polishing step 2
process medaka {

    conda '/home/ubuntu/miniconda3/envs/medaka'

    publishDir "${params.outdir}/06-medaka-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, "${sampleID}.fasta", "${sampleID}.fastq.gz" from medaka

    output:
    set sampleID, "${sampleID}.contigs.racon.medaka.fasta", "${sampleID}.fastq.gz" into unknown
// use r941_min_sup_g507 or r941_min_hac_g507 for model
    """
    medaka_consensus \
    -d ${sampleID}.fasta \
    -i ${sampleID}.fastq.gz \
    -o ${sampleID}_medaka_output \
    -t 10 \
    -m r941_min_sup_g507

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}.contigs.racon.medaka.fasta

    """
}

log.info "If nothing happened, did you inlcude '*.fasta' and '*.fastq.gz' in the --genomes and --trimmedReads options respectively? Also, make sure that the basename of genome and reads is the same!"