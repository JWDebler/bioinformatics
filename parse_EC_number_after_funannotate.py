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



if os.path.isfile("enzyme.dat"):
    enzyme_file = Path("enzyme.dat") 
else:
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
            current_name = ''
        if line.startswith('DE'):
            current_name = current_name + line[4:-2].strip().replace(",","%2C")
        enzyme_dict[current_id] = current_name


with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        elements = []
        for element in line:
            elements.append(element)
        if len(elements) < 8:
            continue
        if "eC_number" in elements[8] or "EC_number" in elements[8]or "ec_number" in elements[8] or "Ec_number" in elements[8] :
            #parse EC number
            search = re.search('[eE][cC]_number=(.+?)$', elements[8])
            ec=search.group(1).split(';')[0]
            #get enzyme name
            # if EC number is incomplete (1, 1.1 or 1.1.1 instead of 1.1.1.1), don't bother looking it up
            if ec.count('.') == 3:
                enzyme_name = enzyme_dict[ec]
            else:
                enzyme_name = ''                
            # Check if enzyme name is blank and remove EC_number tag
            if len(enzyme_name) < 1 :
                element_8 = elements[8].split(";")
                idx = 0
                for i in element_8:
                    idx +=1
                    if i.startswith("product"):
                        idx_product = idx
                    if i.startswith("eC_number") or i.startswith("EC_number") or i.startswith("ec_number") or i.startswith("Ec_number"):
                        idx_EC = idx

                element_8[idx_product-1] = 'product=hypothetical protein'
                #remove EC tag
                del element_8[idx_EC - 1]
                #rebuilt column 8 string
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
            
            # if enzyme name contains 'Transferred entry' (moved) remove EC_number tag and find correct one
            elif "Transferred entry" in enzyme_name:
                #some ECs have been transferred (split) into several new ones. I'll just pick the first one...
                ec = (enzyme_name[19:]).split(" ")[0]
                enzyme_name = enzyme_dict[ec]
                #sometimes entries are transferred to an EC which later gets deleted, this should catch those rare cases
                if enzyme_name == "Deleted entry":
                    element_8 = elements[8].split(";")
                    idx = 0
                    idx_product = 0
                    for i in element_8:
                        idx +=1
                        if i.startswith("product"):
                            idx_product = idx
                            del element_8[idx_product - 1]
                    idx = 0
                    for i in element_8:
                        idx += 1
                        if i.startswith("eC_number") or i.startswith("EC_number") or i.startswith("ec_number") or i.startswith("Ec_number"):
                            idx_EC = idx
                            del element_8[idx_EC -1]
                    
                    
                    element_8_new =';'.join(map(str,element_8))
                    elements[8] = element_8_new
                    print(*elements, sep='\t')
                    element_8_new =';'.join(map(str,element_8))


                else:
                    element_8 = elements[8].split(";")
                    idx = 0
                    idx_product = 0
                    for i in element_8:
                        idx +=1
                        if i.startswith("product"):
                            idx_product = idx
                        if i.startswith("eC_number") or i.startswith("EC_number") or i.startswith("ec_number") or i.startswith("Ec_number"):
                            idx_EC = idx

                    element_8[idx_product-1] = 'product='+enzyme_name
                    element_8[idx_EC-1] = 'EC_number='+ec
                    element_8_new =';'.join(map(str,element_8))
                    elements[8] = element_8_new
                    print(*elements, sep='\t')
                    element_8_new =';'.join(map(str,element_8))

            # if enzyme name contains 'Deleted entry' (deleted) remove EC_number tag
            elif "Deleted entry" in enzyme_name:
                
                ec = (enzyme_name[19:])
                enzyme_name = enzyme_dict[ec]
                element_8 = elements[8].split(";")
                idx = 0
                idx_product = 0
                for i in element_8:
                    idx +=1
                    if i.startswith("product"):
                        idx_product = idx
                        del element_8[idx_product - 1]
                idx = 0
                for i in element_8:
                    idx += 1
                    if i.startswith("eC_number") or i.startswith("EC_number") or i.startswith("ec_number") or i.startswith("Ec_number"):
                        idx_EC = idx
                        del element_8[idx_EC -1]
                
                
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
                element_8_new =';'.join(map(str,element_8))
                
            
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
                    if i.startswith("eC_number") or i.startswith("EC_number") or i.startswith("ec_number") or i.startswith("Ec_number"):
                        idx_EC = idx

                element_8[idx_product-1] = 'product='+enzyme_name
                element_8_new =';'.join(map(str,element_8))
                elements[8] = element_8_new
                print(*elements, sep='\t')
        else:
            #print("++++++++++ no ec")
            print(*elements, sep='\t')
            continue
