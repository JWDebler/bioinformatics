#!/usr/bin/perl -w

###INPUT ARGUMENTS ######################################################
my $inputdir = "/p8/mcc_y95/jdebler/canu_input"; 
my $outdir = "/p8/mcc_y95/jdebler/canu_output";

### CHANGE VALUES AS PER NEED ##########################################
my $username = "johannes.debler"; ## curtin university email. Not implemented at the moment
my $jobname = "canu2";
my $nodes = 1;
my $walltime = "48:00:00";
my $modules = "
#module load canu2
module load java64/1.8.0_202
";

my $cmdline = "canu -p PREFIX -d outdir\/PREFIX\ useGrid=false corThreads=48 genomeSize=45m -nanopore input1";


############## DO NOT EDIT ANYTHING BELOW THIS LINE ############

$inputdir =~ s/\/$//; 
$outdir =~ s/\/$//;

# my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);

my $subdir = "$inputdir\/subscripts"; 
my $cmd = $cmdline;

print "\n	______________________________________________________________________\n
	---A slurm script generator for CANU 2.0 by Fredrick M Mobegi, PhD.---
	______________________________________________________________________\n";
print "\n\n	PLEASE FOLLOW ALL INSTRUCTIONS CLEARLY TO ENSURE YOUR JOBS RUN IN ACCORDANCE TO SERVER FAIR-USE RULES!!!\n";
print "\n	Remember to change emailing domain \<\<\@institution.edu.au\>\> in line 144 of this script if NOT working on PAWSEY HPC.\n	The USERNAME is automalically generated :)\n";
print "\n	Creating output directory: 
	$outdir\n";

if (-e $subdir and -d $subdir) {
    print "	The $subdir already exists. No changes will be made.
	all subscripts from other pipelines will be not be affected\n";
} else {
    mkdir $subdir;
}

if (-e $outdir and -d $outdir) {
    print "	Warning: a folder with the name $outdir elready exists!! Data may be overwritten.\n";
} else {
    mkdir $outdir;
    print "$outdir created successfully!";
}

print "\n	The output directory for assembly files is:
	$outdir\n";

opendir(DIR,$inputdir) or die print "Provide input folder containing interleaved pair-end or single-end fasta files for taxonomic profiling.\n$!"; 
my %files=(); 
my %file_PREFIX=();

while(my $file = readdir(DIR)){
		if ($file =~ /.fna$/ or $file =~ /.fa$/ or  $file =~ /.fasta$/ or  $file =~ /.fastq.gz$/) {
			if ($file =~ m/(\w+)./){
				my $PREFIX = $1; 
				#$PREFIX =~ s/_[L|R].*//;
				$file_PREFIX{$PREFIX}++; 
				$files{$PREFIX}{name}{$file}=1;
			}
		}
	}
my $c=0; my @files_sub=(); 

foreach my $t1(keys%file_PREFIX){ 
	my @file_name=(); 
	$c++; push 
	@file_name,"$inputdir\/$_" foreach (keys%{$files{$t1}{name}});
	my @sample_name = sort @file_name;
#	print $t1,"\t",$sample_name[0],"\t",$sample_name[1],"\n";

    $cmdline =~ s/outdir/$outdir/;
    $cmdline =~ s/PREFIX/$t1/g;
    $cmdline =~ s/input1/$sample_name[0]/;
    $cmdline =~ s/outdir/$outdir/; 
    #$cmdline =~ s/PREFIXb/$t1/; 
    #$cmdline =~ s/PREFIXc/$t1/;

    ################ START SBATCH CODE #################
my $line = "\#!/bin/bash

### Job Name
#SBATCH --job-name=$jobname-$c

### Set email type for job
### Accepted options: NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-type=FAIL

### email address for user
#SBATCH --mail-user=$username\@curtin.edu.au

### Queue name that job is submitted to
#SBATCH --partition=y95

### Request nodes
#SBATCH --nodes=$nodes
#SBATCH --time=$walltime
#SBATCH --export=NONE

#export OMP_NUM_THREADS=24

echo Running on host `hostname`
echo Time is `date`

#module(s) if required module load application_module
$modules
export PATH=\"\$PATH\:\/p8\/mcc_y95\/tools\/canu-2.0\/Linux-amd64/bin\/\"

#Run the executable
$cmdline \n
";

    open (SUBFILE,">$subdir\/$jobname\_$c\.sub"); 
    print SUBFILE $line,"\n";
    push @files_sub,"$subdir\/$jobname\_$c\.sub";
    $cmdline= $cmd;
}
################## END SBATCH CODE ####################


my $num_args = $#ARGV + 1;

if ($num_args !=1){
print "\n	Your sbatch submission scripts are located at\:\n	$inputdir\/subscripts
	Inspect a few of these scripts carefully and make sure they have no errors.
	The scripts should look similar to any \<sbatch.sub\> script you could have generated manually to run a single job.
	+++++++++++++++++++++++++++++++++\n
	To submit all tasks\, re-run this script as follows\:
 
	perl multi_slurmScripts_generator_CANU-2.0.pl \-\-run

	All done!
	+++++++++++++++++++++++++++++++++\n\n";
exit;
} 
elsif ($ARGV[0] eq "--run"){	
	# map { print "$_\n" } @ARGV;
	print "	Your jobs will now be submitted....\n\n	Run squeue -u $username to follow on their progress!\n\n";
	foreach my $n (0 ..scalar(@files_sub)-1){
		my $sub_script = $files_sub[$n];
		`sbatch $sub_script`;
	}
} 
else {
	print "Please provide the correct commands for this script to work!!\n";
	exit;
}

