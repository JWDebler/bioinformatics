import requests
from bs4 import BeautifulSoup

import os
import urllib.request
import csv
import argparse
import re
from pathlib import Path

# This parses the eC_number from the output of funannotate and pulls the product name from a database and adds it as the 'product' tag

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    input_file = Path("test.gff3")
    print("No input file provided, use '-i' and supply a .gff file")
    #raise SystemExit

# def scrape_ec(ec):
#     url = "https://enzyme.expasy.org/EC/" + ec
#     html_content = requests.get(url).text
#     requests.session().close()
#     soup = BeautifulSoup(html_content, "lxml")
#     table = soup.find_all('table')[0]
#     rows = []
#     for child in table.children:
#         row = []
#         for td in child:
#             try:
#                 row.append(td.text.replace('\n', ''))
#             except:
#                 continue
#         if len(row) > 0:
#             rows.append(row)

#     enzyme_name_raw = rows[1][0]
#     enzyme_name = enzyme_name_raw[0:-1].replace(",","%2C")
    
#     # old attempt
#     #ecnumber = soup.title.prettify()
#     #enzyme_name = (ecnumber[(ecnumber.find(ec))+len(ec)+1:-10]).replace(",","%2C")# same change for existing names in GFF
#     print(ec, enzyme_name)
#     return enzyme_name

#download https://ftp.expasy.org/databases/enzyme/enzyme.dat
#and parse into a dictionary
enzyme_url = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
urllib.request.urlretrieve(enzyme_url, "enzyme.dat")

current_id = ''
current_name = ''
enzyme_dict = {}
enzyme_file = Path("enzyme.dat") 

with open(enzyme_file) as file:
    lines = file.readlines()
    for line in lines:
        if line.startswith('ID'):
            current_id = line[4:].strip()
        if line.startswith('DE'):
            current_name = line[4:-2].strip().replace(",","%2C")
        enzyme_dict[current_id] = current_name


with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        elements = []
        for element in line:
            elements.append(element)
        if len(elements) < 8:
            continue
        if "eC_number" in elements[8] or "EC_number" in elements[8] :
            #parse EC number
            search = re.search('[eE]C_number=(.+?)$', elements[8])
            ec=search.group(1).split(';')[0]
           
            #get enzyme name
            # if EC number is incomplete (1, 1.1 or 1.1.1 instead of 1.1.1.1), don't bother looking it up
            if ec.count('.') == 3:
                enzyme_name = enzyme_dict[ec]
                #enzyme_name = scrape_ec(ec) #old method
                #print('scrape =====> all good') #debug
            else:
                enzyme_name = ''                
                #print('=====> incomplete ec') #debug
            # Check if enzyme name is blank and remove EC_number tag
            if len(enzyme_name) < 1 :
                element_8 = elements[8].split(";")
                idx = 0
                for i in element_8:
                    idx +=1
                    if i.startswith("eC_number") or i.startswith("EC_number"):
                        break
                del element_8[idx - 1]
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
            
            # if enzyme name contains 'entry' (moved or deleted) remove EC_number tag
            elif "entry" in enzyme_name:
                element_8 = elements[8].split(";")
                idx = 0
                for i in element_8:
                    idx +=1
                    if i.startswith("eC_number") or i.startswith("EC_number"):
                        break
                del element_8[idx - 1]
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
            
            # otherwise update product name
            else:
                #print(enzyme_name) #debug
                element_8 = elements[8].split(";")
                idx = 0
                idx_product = 0
                for i in element_8:
                    idx +=1
                    if i.startswith("product"):
                        idx_product = idx
                    if i.startswith("eC_number") or i.startswith("EC_number"):
                        idx_EC = idx

                element_8[idx_product-1] = 'product='+enzyme_name
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
        else:
            print(*elements, sep='\t')
            continue
