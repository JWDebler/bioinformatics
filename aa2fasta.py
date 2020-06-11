import os
import csv
import argparse
import glob
from pathlib import Path

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.aa" file')
parser.add_argument('-o', '--output', help='path to an output file (default = inputfilename.fasta)')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
    print(input_file)
else:
    input_file_name = '*.aa'
    input_file_path= os.path.join(os.getcwd(), input_file_name)
    for element in glob.glob(input_file_path):
        input_file = Path(element)

if args.output:
    output_file = Path(args.output)
else:
    output_file = Path(os.path.join(os.getcwd(),os.path.basename(input_file)+'.fasta'))

fasta = {}

print("Reading: ", input_file)

#If file with same outputname exists, delete it first
if os.path.isfile(output_file):
    print("Outputfile already exists, overwriting")
    os.unlink(output_file)


with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        #If something is in the first colum, use that as the identifier 
        if line[0] != "":
            id = line[0]
            fasta[id] = ''
            #Some protein sequences contain the star symbol for stop codon, but interproscan crashes if this is still there
            if line[1][-1] == "*":
                fasta[id] += line[1][:-1]
            else:
                fasta[id] = line[1]
        elif line[1] != "":
            #Some protein sequences contain the star symbol for stop codon, but interproscan crashes if this is still there
            if line[1][-1] == "*":
                fasta[id] += line[1][:-1]
            #For some reason some last lines contain a zero. This needs to go
            elif line[1] != '0':
                fasta[id] += line[1]

print("Writing: ", output_file)
for element in fasta:
    print('>'+element, '\n'+fasta[element], file=open(output_file,"a"))