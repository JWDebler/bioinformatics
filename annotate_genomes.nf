#!/usr/bin/env nextflow

params.genomes = false
params.outdir = "results"
params.genemark = '/opt/genemark-ES/gmes_petap.pl'

if ( params.genomes ) {
    genomes = Channel
    .fromPath(params.genomes, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{genomesFordos2unix}
} else {
    log.info "No genomes supplied."
    exit 1
}



//Making sure there are no strange invisible characters left that would mess with programs down the road
//and also cleaning up fna headers with 'sed', because they contain spaces which mess
//with other tools down the road

process dos2unix {
    tag {sampleID}

    input:
    set sampleID, 'genome.fasta' from genomesFordos2unix

    output:
    set sampleID, "${sampleID}.fasta" into dos2unixed

    """
    dos2unix -n genome.fasta ${sampleID}.fasta
    """
}

process addSpeciesNameTofastaHeadersContigs {
    tag {sampleID}
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set sampleID, "${sampleID}.fasta" from dos2unixed

    output:
    set sampleID, "${sampleID}.fasta" into contigsForGenemark

    """
    sed 's,>,>${sampleID}.,g' -i ${sampleID}.fasta
    sed 's, .*\$,,g' -i ${sampleID}.fast  

    """
}

//GenemarkES annotation
process annotation_genemark {
    tag {sampleID}
    cpus 14
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.gtf'

    input:
    set sampleID, "${sampleID}.fasta" from contigsForGenemark

    output:
    set sampleID, "${sampleID}.genemark.gtf", "${sampleID}.fasta" into genemarkOutput

    """
    ${params.genemark} --ES --fungus --cores ${task.cpus} --sequence ${sampleID}.fasta
    mv genemark.gtf ${sampleID}.genemark.gtf
    """
}

//pull proteins out of genemark annotation
process extractProteinsFromGenemark {
  tag {sampleID}
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.proteins.fasta'

  input:
  set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkOutput

  output:
  set sampleID, "${sampleID}.genemark.proteins.fasta" into proteinsFromGenemark

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.fasta
  """
}

