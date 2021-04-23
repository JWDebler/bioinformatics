def helpMessage() {
    log.info"""
    # Protein annotation using deepSig, EffectorP2.0 and Interproscan (including SignalP4.1)
    A pipeline to annotate proteins with signal peptides, EffectorP Score and via interproscan

    ## Examples
    nextflow run protein_annotation.nf \
    --proteome "proteome/*.fasta" \
    

    ## Parameters
    --proteome <glob>
        Required
        A glob of the fasta proteome to be annotated.
        The basename of the file is used as the proteome name.

    --outdir <path>
        Default: `proteome_annotation`
        The directory to store the results in.

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.proteome = false
params.outdir = "proteome_annotation"

if ( params.proteome ) {
    proteome = Channel
    .fromPath(params.proteome, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { unlabledProteins }
} else {
    log.info "No proteome supplied."
    exit 1
}

//GeneMark doesnt have the species ID in the fasta header, therefor we add it here
process appendSpeciesID {
  tag { id }
  publishDir "${params.outdir}", mode: 'copy'

  input:
  set id, "${id}.fasta" from unlabledProteins

  output:
  set id, "${id}.proteins.fasta" into labledProteins

  """
  awk '/^>/ {sub(/>/, ">${id}.", \$1); print \$1} /^[^>]/ {print \$0}' ${id}.fasta > ${id}.proteins.fasta
  """
}

labledProteins
.tap{proteinsForDeepsig}
.tap{proteinsForInterproscan}
.tap{proteinsForEffectorP}


process effectorP {
  tag { id }
  publishDir "${params.outdir}", mode: 'copy'

  input:
  set id, "${id}.input.fasta" from proteinsForEffectorP

  output:
  file "${id}.effectorP.tsv"

  """
  cp ${id}.input.fasta ${id}.local.input.fasta
  docker run -v \$(pwd):/data/ taniguti/effectorp2 python EffectorP_2.0/Scripts/EffectorP.py \
  -s \
  -o /data/${id}.effectorP.tsv \
  -i /data/${id}.local.input.fasta
  """
}

process deepsig {
  tag { id }
  publishDir "${params.outdir}", mode: 'copy'

  input:
  set id, "${id}.fasta" from proteinsForDeepsig

  output:
  file "${id}.deepsig.out"

  """
  cp ${id}.fasta ./${id}.local.fasta 
  docker run -v \$(pwd):/data/ bolognabiocomp/deepsig -f ${id}.local.fasta  -k euk -o ${id}.deepsig.out
  """

}

process interproscan {
  tag { id }
  publishDir "${params.outdir}", mode: 'copy'
  cpus 12

  input:
  set id, "${id}.proteins.fasta" from proteinsForInterproscan

  output:
  file "${id}.interproscan.tsv"

// currently excluding MobiDBLite, as it fails with python3.8 and TMHMM as there is a path problem somewhere
  """
  /opt/interproscan/current/interproscan.sh \
  --cpu ${task.cpus} \
  --seqtype p \
  --disable-precalc \
  --goterms \
  --pathways \
  --iprlookup\
  --input ${id}.proteins.fasta \
  --output-file-base ${id}.interproscan \
  --format tsv
  """

}
