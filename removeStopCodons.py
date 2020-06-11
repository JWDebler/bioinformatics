import os
import csv
import argparse
import glob
from pathlib import Path

#Removes '*' stop codons from fasta files to prevent downstream tools like Interproscan from crashing.
#If stop codon is detected at end of sequence, it just removes it, if detected in the middle, 
#it removes the entire protein as that is a hint of translation in the wrong reading frame.

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.aa" file')
parser.add_argument('-o', '--output', help='path to an output file (default = inputfilename.fasta)')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
    print(input_file)
else:
    input_file_name = '*.fasta'
    input_file_path= os.path.join(os.getcwd(), input_file_name)
    for element in glob.glob(input_file_path):
        input_file = Path(element)

if args.output:
    output_file = Path(args.output)
else:
    output_file = Path(os.path.join(os.getcwd(),os.path.basename(input_file)+'.fasta'))

fasta = {}
print("=================================================================")
print("Reading: ", input_file)
print("=================================================================")

#If file with same outputname exists, delete it first
if os.path.isfile(output_file):
    print("=================================================================")
    print("Outputfile already exists, overwriting")
    print("=================================================================")
    os.unlink(output_file)

#Parsing fasta file
EndStopCodons = 0
InlineStopCodons = 0

with open(input_file) as file:
    input = file.read().splitlines()
    for line in input:
        if not line.strip(): continue
        if line[0] == '>':
            id = line[1:]
            fasta[id] = ''
        else:
            if line[-1] == '*':
                fasta[id] += line[:-1]
                EndStopCodons += 1
            else:
                fasta[id] += line

delete = [key for key in fasta if '*' in fasta[key]]

for key in delete: 
    del fasta[key]
    InlineStopCodons += 1

print("Removed", EndStopCodons, "stop codons from end of protein.")
print("Removed", InlineStopCodons, "sequence(s) due to inline stop codon(s).")
print("=================================================================")
print("Writing: ", output_file)
print("=================================================================")
for element in fasta:
    print('>'+element, '\n'+fasta[element], file=open(output_file,"a"))