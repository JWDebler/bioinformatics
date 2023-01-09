import os
import argparse
import re
from pathlib import Path

# I had to modify telomeric regions in my genome and needed a quick way to offset all the gene coding features
# by the number of bases I added
# Run by supplying a comma separated list of offsets in bp and they're respective contigs.
# Example: Offset contig 1 by 50 bp and contig 2 by -35 bp:
# gff3_offset.py -o 50,-35 -c ctg01,ctg02 -i input.gff

#parse commandline arguments
parser = argparse.ArgumentParser(description='Parses a GFF file and allows you to offset the location of the elements')
parser.add_argument('-i', '--input', help='path to a GFF file')
parser.add_argument('-o', '--offset', required=True, help="comma separated list of integers, positive for bases added or negative for bases removed")
parser.add_argument('-c', '--contig', required=True, help="comma separated list of contig names that need coordinates offset")

args = parser.parse_args()
#input_file = "test1.gff3"

if args.input:
    input_file = Path(args.input)
else:
    print("No input file provided, use '-i' and supply a GFF file")
    #input_file = "test1.gff3" #testfile
    raise SystemExit

offset = args.offset.split(",")
contig = args.contig.split(",")

#check if same amount of arguments are supplied for contigs and offset
if len(offset) != len(contig):
    print("ERROR, supplied arguments don't match, you supplied",len(offset),"offset(s) and",len(contig),"contig(s).")
    raise SystemExit

# Open the GFF3 file
with open(input_file, "r") as f:
    # Iterate over the lines in the file
    for line in f:
        # Skip lines that start with '#' (these are comments)
        if line[0] == "#" or line[0] == " " or line[0] == "\n":
            continue
        # Split the line on the tab character
        fields = line.strip().split("\t")
        # The first 8 fields are required and contain information about the feature
        seqid, source, type, start, end, score, strand, phase = fields[:8]
        # The last field contains the attributes for the feature
        attributes = fields[8]

        if seqid in contig:
            location = contig.index(seqid)
            #print(seqid,location)
            start = int(start) + int(offset[location])
            end = int(end) + int(offset[location])
            print("\t".join([seqid, source, type, str(start), str(end), score, strand, phase, attributes]))

        else: 
            continue