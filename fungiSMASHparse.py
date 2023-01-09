#!/usr/bin/env python3 

import os
import argparse
import re
from pathlib import Path

# This script extracts the secondary metabolite clusters predicted by fungiSMASH 
# https://fungismash.secondarymetabolites.org/
# It outputs them in BED (default) or GFF3 format
# This is now updated for fungiSMASH 7 index.html as formatting seems to have changed from some spaces to tabs

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
    input_file = "index.html"
    #raise SystemExit

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
    #print(content[element].split(' ')[7])
    chromosome = content[element].split(' ')[2] #P9424_ctg_01
    clustertype = content[element].split(' ')[7] #T1PKS etc

    # find line that has "Location: x - x ..." and remove the commas from the number
    start = content[element+8].split(' ')[1].replace(",","")
    stop = content[element+8].split(' ')[3].replace(",","")

    if outputformat == 1 :
        print(chromosome + "\t" + start + "\t" + stop + "\t" + clustertype)
    else:
        print(chromosome + "\t" + "fungiSHMASH7" + "\t" + "cluster" + "\t" + start + "\t" + stop + "\t" + "." + "\t" + "+" + "\t" + "0" +  "\t" + "NAME="+clustertype)
