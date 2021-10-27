
// before running this script, you need to manually concatenate the demultiplexed fastq files. 
// This script expects one fastq file per genome. guppy can do the barcode demultiplexing during basecalling
// now and creates lots of small fastq (or fastq.gz) files in folders called 'barcodeXX', where
// 'XX' stands for the barcode number i.e. 'barcode06'
// For this script to work and name everything correctly you need to concatenate all those files into one .fastq 
// file, NOT fastq.gz!
// There might be a problem of read duplication. To be sure run this bit of code over the concatenated fastq files
//
// for sample in `ls *.fastq | cut -f1 -d'.'`; do cat $sample.fastq | seqkit rmdup -n -o $sample.clean.fastq; done

def helpMessage() {
    log.info"""
    # Nanopore genome polishing
    A pipeline for polishing genomes assembled from Oxford Nanopore reads using Racon and Medaka.

    ## Examples
    nextflow run nanopore_polishing.nf \
    --nanoporeReads "03-trimmed-fastq/*.fastq.gz"
    

    ## Parameters
    --nanoporeReads <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --pacbioReads <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --illuminaReads <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --outdir <path>
        Default: `assembly`
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

params.cores = 14
params.nanoporeReads = false
params.outdir = "assembly"
params.pilon = "/opt/pilon/pilon-1.24.jar"

if ( params.nanoporeReads ) {
    nanoporeReads = Channel
    .fromPath(params.nanoporeReads, checkIfExists: true, type: "file")
    .map {it -> [it.simpleName, it]}
    .tap { nanoporeReadsForAssembly }
} else {
    log.info "No nanopore reads supplied, did you include '*.fastq.gz'?"
    exit 1
}

if ( params.pacbioReads ) {
    nanoporeReads = Channel
    .fromPath(params.pacbioReads, checkIfExists: true, type: "file")
    .map {it -> [it.simpleName, it]}
    .tap { pacbioReadsForAssembly }
} else {
    log.info "No pacbio reads supplied, did you include '*.fastq.gz'?"
    exit 1
}

if ( params.illuminaReads ) {
    illuminaReads = Channel
    .fromPath(params.illuminaReads, checkIfExists: true, type: "file")
    .map {it -> [it.simpleName, it]}
    .groupTuple(sort: true)
    .map {sid, li -> [sid, li[0], li[1], li[2]]}
    .tap {illuminaReadsForPolishing}
} else {
    log.info "No illumina reads supplied, did you include '*.fastq.gz'?"
    exit 1
}

nanoporeReadsForAssembly
.join(pacbioReadsForAssembly)
.tap { readsForAssembly }


process versions {
    publishDir "${params.outdir}/", mode: 'copy'

    output:
    file 'versions.txt'

    """
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    echo racon: >> versions.txt
    racon --version >> versions.txt
    echo --------------- >> versions.txt
    echo hisat2 >> versions.txt
    hisat2 --version | head -n 1 >> versions.txt
    echo --------------- >> versions.txt
    echo pilon >> versions.txt
    java -jar ${params.pilon} --version >> versions.txt

    """

}

// genome assembly
process Canu {
    tag {sampleID}
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta.gz'
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.report'

    memory '30 GB'

    input:
    tuple sampleID, 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz' from readsForAssembly

    output:
    tuple sampleID, "${sampleID}.contigs.fasta", 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz' into racon
    tuple sampleID, "${sampleID}.correctedReads.fasta.gz" into correctedReads
    file "${sampleID}.canu.report"

    """
    canu \
    -p ${sampleID} \
    -d ${sampleID} \
    genomeSize=45m \
    minInputCoverage=5 \
    stopOnLowCoverage=5 \
    -fast \
    -nanopore input.nanopore.fastq.gz \
    -pacbio input.pacbio.fastq.gz

    cp ${sampleID}/*contigs.fasta ${sampleID}.contigs.fasta
    cp ${sampleID}/*correctedReads.fasta.gz ${sampleID}.correctedReads.fasta.gz
    cp ${sampleID}/*.report ${sampleID}.canu.report
    """

}

// polishing step 1
process racon {
    tag {sampleID}
    publishDir "${params.outdir}/05-racon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    tuple sampleID, 'input.fasta', 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz' from racon

    output:
    tuple sampleID, "${sampleID}.contigs.racon.fasta", 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz' into medaka

    """
    minimap2 \
    input.fasta \
    input.nanopore.fastq.gz > minimap.racon.paf

    racon -m 8 -x -6 -g -8 -w 500 -t 14\
    --no-trimming \
    input.nanopore.fastq.gz \
    minimap.racon.paf \
    input.fasta > ${sampleID}.contigs.racon.fasta

    """
}

// polishing step 2
process medaka {

    conda '/home/ubuntu/miniconda3/envs/medaka'

    tag {sampleID} 
    publishDir "${params.outdir}/06-medaka-polish", mode: 'copy', pattern: '*.fasta'

    input:
    tuple sampleID, 'input.fasta', 'input.fastq.gz' from medaka

    output:
    tuple sampleID, "${sampleID}.contigs.racon.medaka.fasta", 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz' into pilon

    """
    medaka_consensus \
    -d input.fasta \
    -i input.nanopore.fastq.gz \
    -o ${sampleID}_medaka_output \
    -m r941_min_sup_g507

    """
}

process hisat2_and_pilon {
  tag {sampleID}
    publishDir "${params.outdir}/06-pilon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    tuple sampleID, 'input.fasta', 'input.nanopore.fastq.gz', 'input.pacbio.fastq.gz', 'fwd.fastq.gz', 'rev.fastq.gz', 'unpaired.fastq.gz' from pilon.join(illuminaReadsForPolishing)

    output:
    tuple sampleID, "${sampleID}.racon.medaka.pilon.fasta"

    """
    hisat2-build input.fasta ${sampleID}

    hisat2 -x ${sampleID} -1 fwd.fastq.gz -2 rev.fastq.gz --threads ${params.cores} --max-intronlen 2000 | samtools view -b | samtools sort -o ${sampleID}.paired.bam

    samtools index ${sampleID}.paired.bam

    hisat2 -x ${sampleID} -U unpaired.fastq.gz --threads ${params.cores} --max-intronlen 2000 | samtools view -b | samtools sort -o ${sampleID}.unpaired.bam

    samtools index ${sampleID}.unpaired.bam

    java -jar ${params.pilon} --genome input.fasta --frags ${sampleID}.paired.bam --unpaired ${sampleID}.unpaired.bam --threads ${params.cores}

    seqkit sort -lr pilon.fasta > pilon.sorted.fasta
    seqkit replace -p '.+' -r '${sampleID}_ctg_{nr}' --nr-width 2 pilon.sorted.fasta > ${sampleID}.racon.medaka.pilon.fasta

    """  
}