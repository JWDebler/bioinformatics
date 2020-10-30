
// before running this script, you need to manually concatenate the demultiplexed fastq files. 
// This script expects one fastq file per genome. guppy can do the barcode demultiplexing during basecalling
// now and creates lots of small fastq (or fastq.gz) files in folders called 'barcodeXX', where
// 'XX' stands for the barcode number i.e. 'barcode06'
// For this script to work and name everything correctly you need to concatenate all those files into one .fastq 
// file, NOT fastq.gz!
// There might be a problem of read duplication. To be sure run this bit of code over the concatenated fastq files
//
// for sample in `ls *.fastq | cut -f1 -d'.'`; do cat $sample.fastq | seqkit rmdup -n -o $sample.clean.fastq; done

params.nanopore = "/data/2020-10-30_Kewell_P9424_canu211/nanopore/*.fastq.gz"
params.illumina = "/data/2020-10-30_Kewell_P9424_canu211/*.fastq.gz"

nanopore_reads = Channel
    .fromPath(params.nanopre, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}

process versions {
    publishDir "${params.outputdir}/", mode: 'copy'

    output:
    file 'versions.txt'

    """
    echo qcat: >> versions.txt
    qcat --version >> versions.txt
    echo --------------- >> versions.txt
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
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

// genome assembly
process Canu {

    publishDir "${params.outputdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta'

    memory '30 GB'

    input:
    set sampleID, 'input.fastq.gz' from nanopore_reads

    output:
    set sampleID, "${sampleID}.contigs.fasta", 'input.fastq.gz' into racon

    """
    canu \
    -p ${sampleID} \
    -d ${sampleID} \
    genomeSize=45m \
    -nanopore input.fastq.gz

    cp ${sampleID}/*contigs.fasta ${sampleID}.contigs.fasta

    """

}

// polishing step 1
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

    racon -m 8 -x -6 -g -8 -w 500 \
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
    set sampleID, "${sampleID}.contigs.racon.medaka.fasta", 'input.fastq.gz' into pilon

    """
    . /home/ubuntu/medaka/venv/bin/activate
    medaka_consensus \
    -d input.fasta \
    -i input.fastq.gz \
    -o ${sampleID}_medaka_output \
    -m r941_min_high_g360

    cp ${sampleID}_medaka_output/consensus.fasta ${sampleID}.contigs.racon.medaka.fasta

    deactivate
    """
}

