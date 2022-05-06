def helpMessage() {
    log.info"""
    # Gene annotaion with GeneMark-Es
    A pipeline to do de novo gene annotation with GeneMark-ES

    ## Examples
    nextflow run annotate_genomes.nf \
    --genomes "genomes/*.fasta"
    

    ## Parameters
    --genome1 <glob>
        Required
        Reference genome fasta file

    --genome2 <glob>
        Required
        Query genome fasta file


    --cores <int>
        Optional
        Default: 14

    --outdir <path>
        Default: `results_genemark`
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


params.genome1 = false
params.genome2 = false
params.outdir = "plotsr_by_contig"
params.cores = 14

if ( params.genome1 ) {
    genomes = Channel
    .fromPath(params.genome1, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{reference}
} else {
    log.info "No reference supplied."
    exit 1
}

if ( params.genome2 ) {
    genomes = Channel
    .fromPath(params.genome2, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{query}
} else {
    log.info "No query supplied."
    exit 1
}

reference
.merge(query)
.set {reference_and_query}


process extract_contigs {
    tag {sampleID}

    input:
    tuple sample1, 'reference.fasta' , sample2, 'query.fasta' from reference_and_query

    output:
    tuple sample1, "reference.fasta.split/*.fasta" into sample1_split
    tuple sample2, "query.fasta.split/*.fasta" into sample2_split

    """
    seqkit split -i reference.fasta
    seqkit split -i query.fasta
    """
}

s1 = sample1_split.flatMap { k, f -> f.collect { fs -> [k, fs] } }.map {sampleID, filename -> [sampleID, filename.baseName.tokenize('_')[-1], filename]}
s2 = sample2_split.flatMap { k, f -> f.collect { fs -> [k, fs] } }.map {sampleID, filename -> [sampleID, filename.baseName.tokenize('_')[-1], filename]}


s1.merge(s2)
.set {minimap_input}

process minimap {

    tag {sampleID}
    publishDir "${params.outdir}/", mode: 'copy', pattern: '*.bam'

    input:
    tuple sample1, id1, fasta1, sample2, id2, fasta2 from minimap_input

    output:
    tuple sample1, id1, fasta1, sample2, id2, fasta2, "${sample1}_${sample2}_ctg_${id1}.bam" into minimap_output

    """
    minimap2 -ax asm5 -t ${params.cores} --eqx ${fasta1} ${fasta2} | samtools sort -O BAM - > ${sample1}_${sample2}_ctg_${id1}.bam
    """
}

process syri {

    conda '/home/ubuntu/miniconda3/envs/plotsr'
    tag {sampleID}
    publishDir "${params.outdir}/", mode: 'copy', pattern: '*.out'

    input:
    tuple sample1, id1, fasta1, sample2, id2, fasta2, "${sample1}_${sample2}_ctg_${id1}.bam" from minimap_output

    output:
    tuple sample1, id1, fasta1, sample2, id2, fasta2, "${sample1}_${sample2}_ctg_${id1}_syri.out" into syri_output

    """
    /opt/syri/syri/bin/syri -c ${sample1}_${sample2}_ctg_${id1}.bam -r ${fasta1} -q ${fasta2} -F B --prefix ${sample1}_${sample2}_ctg_${id1}_

    """
}

process plotsr {

    conda '/home/ubuntu/miniconda3/envs/plotsr'
    tag {sampleID}
    publishDir "${params.outdir}/", mode: 'copy', pattern: '*.png'

    input:
    tuple sample1, id1, fasta1, sample2, id2, fasta2, "${sample1}_${sample2}_ctg_${id1}_syri.out" from syri_output

    output:
    file "${sample1}_${sample2}_ctg_${id1}.png"

    """
    echo ${fasta1}'\t'${sample1}'\t'lw:1.5> genomes.txt
    echo ${fasta2}'\t'${sample2}'\t'lw:1.5 >> genomes.txt

    plotsr --sr ${sample1}_${sample2}_ctg_${id1}_syri.out -W 20 -H 2 --genomes genomes.txt -o ${sample1}_${sample2}_ctg_${id1}.png

    """


}