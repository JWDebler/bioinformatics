import requests
from bs4 import BeautifulSoup

import os
import csv
import argparse
import re
from pathlib import Path

# Parses a list of pfamIDs and returns the domain family name. I made this to check the pfamIDs list at
# https://github.com/darcyabjones/pante/blob/master/data/pfam_ids.txt for deprecated IDs.

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    
    print("No input file provided, use '-i' and supply a .gff file")
    #input_file = "pfamids.txt"
    raise SystemExit

def scrape_pfamid(pfamid):
    url = "https://pfam.xfam.org/family/" + pfamid
    html_content = requests.get(url).text
    soup = BeautifulSoup(html_content, "lxml")
    pfam_title = soup.title.prettify()
    #parse only the name out of the title
    #print(pfam_title)
    pfam_title_start=pfam_title.find("Family:")
    pfam_title_end=pfam_title.find(pfamid)
    return(pfam_title[pfam_title_start:pfam_title_end+len(pfamid)+1])



with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        pfamid = line[0]
        print(scrape_pfamid(pfamid))