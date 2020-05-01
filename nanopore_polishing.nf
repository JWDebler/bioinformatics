

params.basedir = '/home/ubuntu/2020-04-02_nanopore_basecalling/'
params.sequencefiles = "${params.basedir}03-trimmed-fastq/*.fastq.gz"
params.assemblies = "${params.basedir}04-canu-assembly/*.fasta"


rawnanoporereads = Channel
.fromPath(params.sequencefiles)
.map { file -> [file.getSimpleName(), file]}

assemblies = Channel
.fromPath(params.assemblies)
.map { file -> [file.getSimpleName(), file]}

assemblies
.combine(rawnanoporereads, by: 0)
.set { racon }


process racon {

    publishDir "${params.outputdir}/05-racon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from racon

    output:
    set sampleID, "${sampleID}.contigs.racon.fasta", 'input.fastq.gz' into medaka

    """
    minimap2 \
    input.fasta \
    input.fastq.gz > minimap.racon.paf

    racon \
    input.fastq.gz \
    minimap.racon.paf \
    input.fasta > ${sampleID}.contigs.racon.fasta

    """
}

// polishing step 2
process medaka {

    publishDir "${params.outputdir}/06-medaka-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from medaka

    output:
    set sampleID, "${sampleID}.contigs.racon.medaka.fasta", 'input.fastq.gz' into unknown

    """
    medaka_consensus \
    -d input.fasta \
    -i input.fastq.gz \
    -o ${sampleID}_medaka_output \
    -m r941_min_high_g351

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}.contigs.racon.medaka.fasta

    """
}