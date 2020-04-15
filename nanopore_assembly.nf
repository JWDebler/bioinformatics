params.sequencefiles = '/ppgdata/johannes/2020-04-02-nanopore-lentis/FastQ/*/*.fastq'
params.outputdir =     '/home/johannes/rdrive/Ascochyta_genomics-KAMPHL-SE07477/johannes/notebook/2020-04-01_first_minION_run/'
params.canu = '/opt/canu/current/Linux-amd64/bin/canu'

// before running this script, you need to manually concatenate the demultiplexed fastq files. 
// This script expects one fastq file per genome. guppy can do the barcode demultiplexing during basecalling
// now and creates lots of small fastq (or fastq.gz) files in folders called 'barcodeXX', where
// 'XX' stands for the barcode number i.e. 'barcode06'
// For this script to work and name everything correctly you need to concatenate all those files into one .fastq 
// file, NOT fastq.gz!

rawnanoporereads = Channel
.fromPath(params.sequencefiles)
.map { file -> [file.getSimpleName(),file]}

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

// barcode and adapter trimming
process qcat_trim {

    publishDir "${params.outputdir}/03-trimmed-fastq", mode: 'copy', pattern: '*.fastq.gz'
    
    input:
    set sampleID, 'input.fastq' from rawnanoporereads

    output:
    set sampleID, "${sampleID}.trimmed.fastq.gz" into canuAssembly

    """
    qcat -f input.fastq -o "${sampleID}.trimmed.fastq" --detect-middle --trim
    gzip -9 ${sampleID}.trimmed.fastq
    """

}

// genome assembly
process Canu {

    publishDir "${params.outputdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta'

    memory '30 GB'

    input:
    set sampleID, 'input.fastq.gz' from canuAssembly

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
