import os
import argparse
import re
from pathlib import Path

# I downloaded a gff file from JGI which only had CDS and exon features, no gene or mRNA. 
# Importing that into geneious displayed all the CDSs as individual entities instead as 
# connected to each other. This script fixes that by adding the appropriate Parent tags and 
# gene features.

#parse commandline arguments
parser = argparse.ArgumentParser(description='Parses a GFF file downloaded from JGI that only has CDS but not gene features')
parser.add_argument('-i', '--input', help='path to a GFF file')

args = parser.parse_args()
#input_file = "jgi.gff"

if args.input:
    input_file = Path(args.input)
else:
    print("No input file provided, use '-i' and supply a GFF file")
    #input_file = "test1.gff3" #testfile
    raise SystemExit

transcripts = {}
transcript_id = 0
# Open the GFF3 file
with open(input_file, "r") as f:
    # Iterate over the lines in the file
    for line in f:
        # Split the line on the tab character
        fields = line.strip().split("\t")
        # The first 8 fields are required and contain information about the feature
        seqid, source, type, start, end, score, strand, phase = fields[:8]
        # The last field contains the attributes for the feature
        attributes = fields[8]
        transcript = attributes.strip().split(";")[0].replace(" ","_").replace("#","_")
        if type == 'CDS':
            exon_start = int(start)
            exon_end = int(end)
            if transcript not in transcripts:
                transcript_id += 1
                transcripts[transcript] = [seqid,exon_start,exon_end,strand]
                print(seqid + "\t" + source + "\t" +  type + "\t" +  start + "\t" +  end + "\t" +  score + "\t" +  strand + "\t" +  phase  + "\t" + "ID="+transcript+";Parent="+transcript+".gene")
            else:
                if exon_start < transcripts[transcript][1]:
                    transcripts[transcript][1] = exon_start
                    print(seqid + "\t" + source + "\t" +  type + "\t" +  start + "\t" +  end + "\t" +  score + "\t" +  strand + "\t" +  phase  + "\t" + "ID="+transcript+";Parent="+transcript+".gene")
                elif exon_end > transcripts[transcript][2]:
                    transcripts[transcript][2] = exon_end
                    print(seqid + "\t" + source + "\t" +  type + "\t" +  start + "\t" +  end + "\t" +  score + "\t" +  strand + "\t" +  phase  + "\t" + "ID="+transcript+";Parent="+transcript+".gene")

    for key,value in transcripts.items():
        print(value[0] + "\t" + "JGI" + "\t" + "gene" + "\t" + str(value[1]) +"\t" + str(value[2]) + "\t" + "." + "\t" + value[3] + "\t" + "0" + "\t" + "ID=" + key + ".gene")
    