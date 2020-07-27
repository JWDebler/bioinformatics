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

    --SRRids <glob>
        Required
        Textfile with SRR IDs, one per line.

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

if ( params.srrids ) {
    srrids = Channel
    .fromPath(params.srrids, checkIfExists: true, type: "file")
    .splitText()
    .map{it -> it.trim()}
    .set{srrIDsTrimmed}
} else {
    log.info "No SRR IDs text file supplied."
    exit 1
}

//Sometimes SRR datasets are for some reason merged (forward and backward read)
//-split-3 separates those reads back
process dumpfastq {
  publishDir "${params.outdir}/01-fastq", mode: 'copy'

  input:
    val id from srrIDsTrimmed

  output:
    //set id, "*.fastq" into fastqDumpForAlignment
    set id, "*_1.fastq", "*_2.fastq" into fastqDumpForAlignment

  """
  fastq-dump -split-3 $id
  """
}

fastqDumpForAlignment
.tap {fastqForFastQC}
.collect()
.flatten()
.collate(3)
//.collate(2)
.set{fastqDumpForAlignmentAll}

process FastQC {
  publishDir "${params.outdir}/02-fastQC", mode: 'copy'

  input:
  set idFastq, "${idFastq}_1.fastq", "${idFastq}_2.fastq" from fastqForFastQC
  //set idFastq, "${idFastq}.fastq" from fastqForFastQC

  output:
  //set idFastq, "${idFastq}_fastqc.html", "${idFastq}_fastqc.zip"
  set idFastq, "${idFastq}_1_fastqc.html", "${idFastq}_1_fastqc.zip", "${idFastq}_2_fastqc.html", "${idFastq}_2_fastqc.zip"

  """
  /opt/FastQC/fastqc  "${idFastq}_1.fastq" "${idFastq}_2.fastq"
  """
}

process indexAssemblyHisat2 {
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
  tag { "${idAssembly} ${idFastq}" }

  input:
    set idAssembly, "genome.fasta", "${idAssembly}.*.ht2", idFastq, "${idFastq}_1.fastq", "${idFastq}_2.fastq" from indexedAssembly.combine(fastqDumpForAlignmentAll)
    //set idAssembly, "genome.fasta", "${idAssembly}.*.ht2", idFastq, "${idFastq}.fastq" from indexedAssembly.combine(fastqDumpForAlignmentAll)

  output:
    set  idAssembly, "genome.fasta", idFastq, "${idAssembly}.${idFastq}.bam" into bamsToFilter

// either use -U, or -1/-2

// -1 ${idFastq}_1.fastq -2 ${idFastq}_2.fastq \
// -U ${idFastq}.fastq \
  """
  hisat2 -x ${idAssembly} \
  -1 ${idFastq}_1.fastq -2 ${idFastq}_2.fastq \
  --threads 10 \
  --max-intronlen 2000 \
  | samtools view -b \
  | samtools sort -n \
  -o ${idAssembly}.${idFastq}.bam

  """
}

process filterBams {
  publishDir "${params.outdir}/04-filteredBams", mode: 'copy', pattern: '*.bam'

  input: 
  set idAssembly, "genome.fasta", idFastq, "aligned.bam" from bamsToFilter 

  output: 
  set idAssembly, "genome.fasta", idFastq, "${idAssembly}.${idFastq}.filtered.sorted.bam" into alignedAndFiltered

  """
  cp aligned.bam aligned.input.bam
  docker run -v \$PWD:/xxx augustus /opt/augustus-3.3.3/bin/filterBam \
  --uniq \
  --in /xxx/aligned.input.bam \
  --out /xxx/filtered.bam
  samtools sort filtered.bam > ${idAssembly}.${idFastq}.filtered.sorted.bam
  """
}

alignedAndFiltered
.tap {bamToHintsInput}
.tap {indexBamInput}



process indexBam {
  publishDir "${params.outdir}/04-filteredBams/", mode: 'copy', pattern: '*.bai'
  tag { "${idAssembly} ${idFastq}" }

  input:
  set  idAssembly, "genome.fasta", idFastq, "assembly.bam" from indexBamInput

  output:
  set  idAssembly, "genome.fasta", idFastq, "${idAssembly}.${idFastq}.filtered.sorted.bam.bai" into indexedBam

  """
  samtools index assembly.bam ${idAssembly}.${idFastq}.filtered.sorted.bam.bai
  """
}


process bamToHints {
  
  tag { "${idAssembly} ${idFastq}" }

  input:
  set idAssembly, "genome.fasta", idFastq, "assembly.bam" from bamToHintsInput

  output:
  set idAssembly, "genome.fasta", idFastq, "${idAssembly}.${idFastq}.gff" into bamToHintsOutput

  """
  cp assembly.bam assembly.input.bam

  docker run -v \$PWD:/xxx augustus /opt/augustus-3.3.3/bin/bam2hints \
  --intronsonly \
  --maxgaplen=10 \
  --minintronlen=15 \
  --maxintronlen=500 \
  --in=/xxx/assembly.input.bam \
  --out=/xxx/${idAssembly}.${idFastq}.gff
  """
}

process findStrand {
   publishDir "${params.outdir}/05-hints/", mode: 'copy', pattern: '*.gff'
  tag { "${idAssembly} ${idFastq}" }

  input:
  set idAssembly, "genome.fasta", idFastq, "introns.gff" from bamToHintsOutput

  output:
  set idAssembly, "genome.fasta", idFastq, "${idAssembly}.${idFastq}.introns.gff" into findStrandOutput

  """
  cp genome.fasta genome.input.fasta
  cp introns.gff introns.input.gff

  docker run -v \$PWD:/xxx augustus perl /root/augustus/docs/tutorial2018/BRAKER_v2.0.4+/filterIntronsFindStrand.pl \
  /xxx/genome.input.fasta \
  /xxx/introns.input.gff \
  --score > "${idAssembly}.${idFastq}.introns.gff"
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
  set idAssembly, genomes, idFastq, "hints*.gff", "genome.fasta" from mergeHints.combine(assembliesForAugustus, by:0)

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
  --cores=14 \
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
  set idAssembly, "bonafide.f.gtf", "genome.fasta" , "bonafide.gb"from extractAA

  output:
  set idAssembly, "${idAssembly}.trainingset.aa" , "bonafide.gb"into trainingsetAA

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
  set idAssembly, "${idAssembly}.trainingset.nonred.aa" , "bonafide.gb"into trainingsetAAnonredundant

  """
  /opt/Augustus/scripts/aa2nonred.pl \
  trainingset.aa \
  ${idAssembly}.trainingset.nonred.aa \
  --cores=14
  """
}

process cleanupAA {

  input:
  set idAssembly, "nonred.aa" , "bonafide.gb"from trainingsetAAnonredundant

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


return


process augustusGTF {
  publishDir "${params.outdir}/$id/augustus", mode: 'copy'
  tag { id }

  input:
  set id, "hints.gff", "input.fasta" from mergedHintsForAugustus


  output:
  set id, "${id}.augustus.gff", "input.fasta" into augustusGTFoutput

  """
docker run augustus augustus \
 --progress=true \
 --gff3=off\
 --softmasking=1 \
 --uniqueGeneId=true \
 --noInFrameStop=true \
 --hintsfile=hints.gff \
 --/augustus/verbosity=4 \
 --species=PnodSN15 \
 --extrinsicCfgFile=${params.config} \
 input.fasta > ${id}.augustus.gff
  """
}

augustusGTFoutput
.tap{augustusToFasta}
.tap{augustusToBed}

process extractFasta {
  publishDir "${params.outdir}/$id/augustus", mode: 'copy'
  tag { id }

  input:
  set id, "${id}.gff", "input.fasta" from augustusToFasta

  output:
  set id, "${id}.proteins.fasta" into augustusFastas

  """
  /opt/augustus/current/scripts/getAnnoFasta.pl --seqfile input.fasta ${id}.gff
  mv ${id}.aa ${id}.proteins.fasta
  mv ${id}.codingseq ${id}.cds.fasta
  """

}