#!/usr/bin/env python3 

import os
import argparse
import re
from pathlib import Path

# This script extracts the secondary metabolite clusters predicted by fungiSMASH 
# https://fungismash.secondarymetabolites.org/
# It outputs them in BED (default) or GFF3 format

#parse commandline arguments
parser = argparse.ArgumentParser(description='Parses fungiSHMASH index.html file for location of clusters. By default looks for index.html file in current directory')
parser.add_argument('-i', '--input', help='path to "index.html" file downloaded from fungiSMASH')
parser.add_argument('-g', '--gff',action='store_true', help='output in GFF format')
parser.add_argument('-b', '--bed',action='store_true', help='output in BED format (default)')

args = parser.parse_args()

input_file = "index.html" #default filename

if args.input:
    input_file = Path(args.input)
else:
    print("No input file provided, use '-i' and supply the index.html file provided by fungiSMASH")
    #input_file = "index.html"
    raise SystemExit

outputformat = 1 # 1 = BED (default)
                 # 2 = GFF

if args.gff:
    outputformat = 2

regions = []
content = [] # paresed HTML file
linenumber = 0

with open(input_file) as file:
    lines = file.readlines()
    for line in lines:
        linenumber += 1
        content.append(line)
        # find region and type of cluster
        if '<div class ="description-container">' in line:
            regions.append(linenumber)

for element in regions:
    chromosome = content[element].split(' ')[8]
    clustertype = content[element].split(' ')[13]
    start = content[element+8].split(' ')[9].replace(",","")
    stop = content[element+8].split(' ')[11].replace(",","")
    if outputformat == 1 :
        print(chromosome,"\t", start, "\t",stop, "\t",clustertype)
    else:
        print(chromosome,"fungiSHMASH","cluster",start,stop,".","+","0","NAME="+clustertype)
