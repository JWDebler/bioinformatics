// Build your own Augustus model from RNAseq data and a reference genome assembly
// according to: "Hoff KJ, Stanke M. Predicting Genes in Single Genomes with AUGUSTUS. Curr Protoc Bioinforma. 2019;65(1):1-54. doi:10.1002/cpbi.57"

// Prerequisits: 
// - Genome assembly
// - RNAseq data

def helpMessage() {
    log.info"""
    # A pipeline to build your own model for Augustus
    based on RNAseq data and a genome assembly

    ## Examples
    nextflow run build_augustus_model_from_RNAseq.nf \
    --genome "genome/*.fasta" \
    --rnaseq ""
    

    ## Parameters
    --genome <glob>
        Required
        A glob of the fasta genome to be used for model building.
        The basename of the file is used as the genome name.

    --cores <int>
        Optional
        Default: 14

    --reads <glob>
        Required
        fastq.gz file of forward reads.

    --outdir <path>
        Default: `results_augustus_model`
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

params.genome = false
params.rnaseq = false
params.srrids = false
params.outdir = "results"
params.cores = 14

if ( params.genome ) {
    genome = Channel
    .fromPath(params.genome, checkIfExists: true, type: "file")
    .map{file -> [file.simpleName, file]}
    .tap{assemblies}
    .tap{assembliesForAugustus}
} else {
    log.info "No genome supplied, did you include '*.fasta'?"
    exit 1
}

if ( params.reads ) {
    fwd_reads = Channel
    .fromPath(params.forward, checkIfExists: true, type: "file")
    .set{forward_reads}
} else {
    params.reads = "/home/johannes/rdrive/Ascochyta_genomics-KAMPHL-SE07477/johannes/notebook/2020-10-08_Rabiei_RNAseq/02-trimmed_reads/*_{1,2}P.fq.gz"
    //log.info "No reads supplied."
    //exit 1
}

reads = 
Channel.fromFilePairs(params.reads)
.map {sampleID, fwdrevreads -> [sampleID.tokenize('_')[0], fwdrevreads]}
.groupTuple()
.map {sampleID, ary -> [sampleID, ary.transpose()]}
.map {sampleID, ary -> [sampleID, ary[0], ary[1]]}
.set{rawReads}

process combine_reads {
  tag {readID}

  publishDir "${params.outdir}/02-combined-reads", mode: 'copy'

  input: 
  set readID, "${readID}.fwd.*.fastq.gz", "${readID}.rev.*.fastq.gz" from rawReads

  output:
        set readID, "${readID}.fwd.fastq.gz", "${readID}.rev.fastq.gz" into fwdrevreads

        """
        cat ${readID}.fwd.*.fastq.gz > ${readID}.fwd.fastq.gz
        cat ${readID}.rev.*.fastq.gz > ${readID}.rev.fastq.gz \
        """
}

process indexAssembly2 {
  tag { id }

  input:
    set id, "genome.fasta" from assemblies

  output:
    set id, "genome.fasta", "*.ht2" into indexedAssembly

  """
  hisat2-build genome.fasta $id
  """
}

process alignToAssemblyhisat2 {
  publishDir "${params.outdir}/03-bams", mode: 'copy', pattern: '*.bam'
  tag { "${idAssembly} ${readID}" }
  memory '30 GB'

  input:
    set idAssembly, "genome.fasta", "${idAssembly}.*.ht2", readID, "fwd.fastq.gz", "rev.fastq.gz" from indexedAssembly.combine(fwdrevreads)

  output:
    set  idAssembly, "genome.fasta", readID, "${idAssembly}.${readID}.bam" into bamsToFilter

// either use -U, or -1/-2
// -1 ${idFastq}_1.fastq -2 ${idFastq}_2.fastq \
// -U ${idFastq}.fastq \
  """
  hisat2 -x ${idAssembly} \
  -1 fwd.fastq.gz -2 rev.fastq.gz \
  --threads ${params.cores} \
  --max-intronlen 2000 \
  | samtools view -b \
  | samtools sort -n \
  -o ${idAssembly}.${readID}.bam

  """
}

process filterBams {
  publishDir "${params.outdir}/04-filteredBams", mode: 'copy', pattern: '*.bam'
  tag { "${idAssembly} ${readID}" }

  input: 
  set idAssembly, "genome.fasta", readID, "aligned.bam" from bamsToFilter 

  output: 
  set idAssembly, "genome.fasta", readID, "${idAssembly}.${readID}.filtered.sorted.bam" into alignedAndFiltered

  """
  cp aligned.bam aligned.input.bam
  docker run -v \$PWD:/xxx augustus /opt/augustus-3.3.3/bin/filterBam \
  --uniq \
  --in /xxx/aligned.input.bam \
  --out /xxx/filtered.bam
  samtools sort filtered.bam > ${idAssembly}.${readID}.filtered.sorted.bam
  """
}

alignedAndFiltered
.tap {bamToHintsInput}
.tap {indexBamInput}



process indexBam {
  publishDir "${params.outdir}/04-filteredBams/", mode: 'copy', pattern: '*.bai'
  tag { "${idAssembly} ${readID}" }

  input:
  set  idAssembly, "genome.fasta", readID, "assembly.bam" from indexBamInput

  output:
  set  idAssembly, "genome.fasta", readID, "${idAssembly}.${readID}.filtered.sorted.bam.bai" into indexedBam

  """
  samtools index assembly.bam ${idAssembly}.${readID}.filtered.sorted.bam.bai
  """
}


process bamToHints {
  
  tag { "${idAssembly} ${readID}" }

  input:
  set idAssembly, "genome.fasta", readID, "assembly.bam" from bamToHintsInput

  output:
  set idAssembly, "genome.fasta", readID, "${idAssembly}.${readID}.gff" into bamToHintsOutput

  """
  cp assembly.bam assembly.input.bam

  docker run -v \$PWD:/xxx augustus /opt/augustus-3.3.3/bin/bam2hints \
  --intronsonly \
  --maxgaplen=10 \
  --minintronlen=15 \
  --maxintronlen=500 \
  --in=/xxx/assembly.input.bam \
  --out=/xxx/${idAssembly}.${readID}.gff
  """
}

process findStrand {
   publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gff'
  tag { "${idAssembly} ${readID}" }

  input:
  set idAssembly, "genome.fasta", readID, "introns.gff" from bamToHintsOutput

  output:
  set idAssembly, "genome.fasta", readID, "${idAssembly}.${readID}.introns.gff" into findStrandOutput

  """
  cp genome.fasta genome.input.fasta
  cp introns.gff introns.input.gff

  docker run -v \$PWD:/xxx augustus perl /root/augustus/docs/tutorial2018/BRAKER_v2.0.4+/filterIntronsFindStrand.pl \
  /xxx/genome.input.fasta \
  /xxx/introns.input.gff \
  --score > "${idAssembly}.${readID}.introns.gff"
  """
}

findStrandOutput
.collect()
.flatten()
.collate(4)
.groupTuple()
.set {mergeHints}

process mergeHints {
  publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gff'
  input :
  set idAssembly, genomes, readID, "hints*.gff", "genome.fasta" from mergeHints.combine(assembliesForAugustus, by:0)

  output:
  set idAssembly, "merged_introns.gff" , "genome.fasta" into genemark

  """
  cat hints* > merged_introns.gff
  """
}

process genemark {
  publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gtf'
  input:
  set idAssembly, "introns.gff", "genome.fasta" from genemark

  output:
  set idAssembly, "introns.gff",  "genome.fasta", "${idAssembly}.gtf" into filterGenemark

  """
  /opt/genemark-ES/gmes_petap.pl \
  --verbose \
  --sequence=genome.fasta \
  --ET=introns.gff \
  --cores=${params.cores} \
  --soft_mask 1000

  mv genemark.gtf ${idAssembly}.gtf

  """
}

process filterGenemark {
  publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gtf'
  input:
  set idAssembly, "introns.gff", "genome.fasta", "genemark.gtf" from filterGenemark

  output:
  set idAssembly, "genome.fasta", "${idAssembly}.bonafide.gtf", "genemark.gtf" into computeFlankingRegions

  """
  cp genemark.gtf genemark.local.gtf
  cp introns.gff introns.local.gff

  docker run -v \$PWD:/xxx augustus perl /root/augustus/docs/tutorial2018/BRAKER_v2.0.4+/filterGenemark.pl \
  --genemark=/xxx/genemark.local.gtf \
  --introns=/xxx/introns.local.gff

  mv *.good.gtf ${idAssembly}.bonafide.gtf
  """
}

process computeFlankingRegions {

  input:
  set idAssembly, "genome.fasta", "bonafide.gtf", "genemark.gtf" from computeFlankingRegions

  output:
  set idAssembly, "genome.fasta",  "${idAssembly}.bonafide.gtf", "tmp.gb" into filterGenesIn_mRNAname

  """
  cp bonafide.gtf bonafide.local.gtf
  cp genemark.gtf genemark.local.gtf
  cp genome.fasta genome.local.fasta

  flankvalue=\$(docker run -v \$PWD:/xxx augustus /root/augustus/scripts/computeFlankingRegion.pl \
  /xxx/bonafide.local.gtf | \
  awk '/^The flanking_DNA value is:/ {print \$5}')

  docker run -v \$PWD:/xxx augustus /root/augustus/scripts/gff2gbSmallDNA.pl \
  /xxx/genemark.local.gtf \
  /xxx/genome.local.fasta \
  \$flankvalue /xxx/tmp.gb

  mv bonafide.gtf ${idAssembly}.bonafide.gtf
  """

}


process filterGenesIn_mRNAname {
    publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gb'

    input:
    set idAssembly, "genome.fasta", "bonafide.gtf", "tmp.gb" from filterGenesIn_mRNAname

    output:
    set idAssembly, "${idAssembly}.bonafide.gb", "bonafide.gtf", "genome.fasta" into extractTranscriptNames

    """
    cp tmp.gb tmp.local.gb
    cp bonafide.gtf bonafide.local.gtf

    docker run -v \$PWD:/xxx augustus /root/augustus/scripts/filterGenesIn_mRNAname.pl \
    /xxx/bonafide.local.gtf \
    /xxx/tmp.local.gb \
    > bonafide.gb

    mv bonafide.gb ${idAssembly}.bonafide.gb
    """
}

process extractTranscriptNames {

  input:
  set idAssembly, "bonafide.gb", "bonafide.gtf", "genome.fasta" from extractTranscriptNames

  output:
  set idAssembly, "bonafide.f.gtf", "genome.fasta", "bonafide.gb" into extractAA

  """
  cat bonafide.gb | awk '/\\/gene/ {g=gensub(/.*gene="([^[:space:]]+)"/, "\\"\\\\1\\"", "g", \$0); print g}' | sort -u > traingenes.lst
  grep -f traingenes.lst -F bonafide.gtf > bonafide.f.gtf
  """
}

process extractAA {

  publishDir "${params.outdir}/06-training/", mode: 'copy', pattern: '*.aa'

  input:
  set idAssembly, "bonafide.f.gtf", "genome.fasta" , "bonafide.gb" from extractAA

  output:
  set idAssembly, "${idAssembly}.trainingset.aa" , "bonafide.gb" into trainingsetAA

  """
  cp bonafide.f.gtf bonafide.local.gtf
  cp genome.fasta genome.local.fasta

  docker run -v \$PWD:/xxx augustus /root/augustus/scripts/gtf2aa.pl \
  /xxx/genome.local.fasta \
  /xxx/bonafide.local.gtf \
  /xxx/${idAssembly}.trainingset.aa
  """
}

process blastAllvsAll {

  publishDir "${params.outdir}/06-training/", mode: 'copy', pattern: '*.aa'
  
  input:
  set idAssembly, "trainingset.aa" , "bonafide.gb"from trainingsetAA

  output:
  set idAssembly, "${idAssembly}.trainingset.nonred.aa" , "bonafide.gb" into trainingsetAAnonredundant

  """
  /opt/Augustus/scripts/aa2nonred.pl \
  trainingset.aa \
  ${idAssembly}.trainingset.nonred.aa \
  --cores=${params.cores}
  """
}

process cleanupAA {

  publishDir "${params.outdir}/06-training/", mode: 'copy', pattern: '*.gb'

  input:
  set idAssembly, "nonred.aa" , "bonafide.gb" from trainingsetAAnonredundant

  output:
  file "${idAssembly}.bonafide.f.gb"

  """
  grep ">" nonred.aa | perl -pe 's/>//' > nonred.lst

  cat bonafide.gb | \

  perl -ne 'if (\$_ =~ m/LOCUS\\s+(\\S+)\\s/) { \$txLocus = \$1;} elsif (\$_ =~ m/\\/gene=\\"(\\S+)\\"/) {\$txInGb3{\$1} = \$txLocus} if(eof()) {foreach (keys %txInGb3) { print "\$_\\t\$txInGb3{\$_}\\n";}}' > loci.lst
  
  grep -f nonred.lst loci.lst | cut -f2 > nonred.loci.lst

  cp bonafide.gb bonafide.local.gb
  docker run -v \$PWD:/xxx augustus /root/augustus/scripts/filterGenesIn.pl \
  /xxx/nonred.loci.lst \
  /xxx/bonafide.local.gb > ${idAssembly}.bonafide.f.gb

  
  """
}
