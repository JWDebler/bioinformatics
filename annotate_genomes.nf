#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    # Gene annotaion with GeneMark-Es
    A pipeline to do de novo gene annotation with GeneMark-ES

    ## Examples
    nextflow run annotate_genomes.nf \
    --genomes "genomes/*.fasta"
    

    ## Parameters
    --genomes <glob>
        Required
        A glob of the fasta genomes to be annotated.
        The basename of the file is used as the genome name.

    --cores <int>
        Optional
        Default: 14

    --outdir <path>
        Default: `results_genemark`
        The directory to store the results in.

    --species <string>
        Default "Alentis"
        Species to use for Augustus

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false
params.outdir = "genome_annotation"
params.genemark = '/opt/genemark-ES/gmes_petap.pl'
params.cores = 14
params.species = "Alentis"

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
    set sampleID, "${sampleID}.fasta" into contigsForGenemark, contigsForAugustus

    """
    sed 's,>,>${sampleID}.,g' -i ${sampleID}.fasta
    sed 's, .*\$,,g' -i ${sampleID}.fasta  

    """
}

//Augustus annotation
process annotation_augustus {
    tag {sampleID}
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.gff3'

    input:
    set sampleID, "${sampleID}.fasta" from contigsForAugustus

    output:
    set sampleID, "${sampleID}.augustus.gff3", "${sampleID}.fasta" into augustusOutput

    """
    cp ${sampleID}.fasta input.fasta

    docker run \
    -v \$PWD:/xxx \
    -v /opt/Augustus/config:/root/augustus/config \
    augustus \
    augustus \
    --species=${params.species} \
    --progress=true \
    --gff3=on \
    --softmasking=1 \
    --uniqueGeneId=true \
    --noInFrameStop=true \
    --/augustus/verbosity=4 \
    /xxx/input.fasta > ${sampleID}.augustus.gff3
    """
}

//GenemarkES annotation
process annotation_genemark {
    tag {sampleID}
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.gtf'

    input:
    set sampleID, "${sampleID}.fasta" from contigsForGenemark

    output:
    set sampleID, "${sampleID}.genemark.gtf", "${sampleID}.fasta" into genemarkProt, genemarkGFF

    """
    ${params.genemark} --ES --fungus --cores ${params.cores} --sequence ${sampleID}.fasta
    mv genemark.gtf ${sampleID}.genemark.gtf
    """
}

//Convert genemark GTF to GFF with fixed CDS phase

process genemarkGTFtoGFF {
    tag {sampleID}
    publishdir "${params.outdir}", mode: 'copy', pattern: '*.gff3'

    input:
    set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkGFF

    output:
    set sampleID, "${sampleID}.genemark.gff3"

    """
    agat_convert_sp_gxf2gxf.pl -g ${sampleID}.genemark.gtf -gvo 3 -o genemark.gff
    agat_sp_manage_IDs.pl -gff genemark.gff --collective -o genemark.CDS.gff
    agat_sp_fix_cds_phases.pl -g genemark.CDS.gff -f input.fasta -o ${sampleID}.genemark.gff3
    """
}

//pull proteins out of genemark annotation
process extractProteinsFromGenemark {
  tag {sampleID}
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.proteins.fasta'

  input:
  set sampleID, "${sampleID}.genemark.gtf", "input.fasta" from genemarkProt

  output:
  set sampleID, "${sampleID}.genemark.proteins.fasta" 

  """
  /opt/genemark-ES/get_sequence_from_GTF.pl ${sampleID}.genemark.gtf input.fasta
  mv prot_seq.faa ${sampleID}.genemark.proteins.fasta
  """
}

//pull proteins out of augustus annotation
process extractProteinsFromAugustus {
  tag {sampleID}
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.proteins.fasta'

  input:
  set sampleID, "input.gff3", "input.fasta" from augustusOutput

  output:
  set sampleID, "${sampleID}.augustus.proteins.fasta" 

  """
  /agat_sp_extract_sequences.pl -g input.gff3 -f input.fasta -p -o ${sampleID}.augustus.proteins.fasta
  """
}