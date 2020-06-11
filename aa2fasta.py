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

with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        if line[0] != "":
            id = line[0]
            fasta[id] = ''
            if line[1][-1] == "*":
                fasta[id] += line[1][:-1]
            else:
                fasta[id] = line[1]
        elif line[1] != "":
            if line[1][-1] == "*":
                fasta[id] += line[1][:-1]
            elif line[1] != '0':
                fasta[id] += line[1]

for element in fasta:
    print('>'+element, '\n'+fasta[element], file=open(output_file,"a"))