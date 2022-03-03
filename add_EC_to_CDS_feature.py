import csv
import argparse
import re
from pathlib import Path

# parses gff file created by 'parse_EC_number_after_funannotate.py' and if found adds 'product' tags to 'CDS'
# as that seems to be what NCBI expects

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    input_file = Path("test1.gff3")
    print("No input file provided, use '-i' and supply a .gff file")
    #raise SystemExit

with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    id_product_mapper = {}
    for line in input:
        if line[0].startswith('#'):
            continue
        elements = []
        
        for element in line:
            elements.append(element)
        
        # fill dictionary with EC number and matching ID
        if "product=" in elements[8]:
            #print(elements[8]) #debug
            product_string = re.search('product=(.+?)$', elements[8])   
            id_string = re.search('ID=(.+?)$', elements[8]) 
            product=product_string.group(1).split(';')[0]
            id=id_string.group(1).split(';')[0]
            #print(product, id) debug
            id_product_mapper[id] = 'product='+product+';'

#print(id_product_mapper)

with open(input_file) as file:   
    input = csv.reader(file, delimiter='\t')
    for line in input:
        if line[0].startswith('#'):
            continue
        elements = []
        for element in line:
            elements.append(element)
        if elements[2] == 'CDS':
            #print('before:', line)     #debug
            id_string = re.search('ID=(.+?)$', elements[8]) 
            id=id_string.group(1).split(';')[0]
            id = id.split('.')[0]
            if id in id_product_mapper:
                elements[8] += ';' + id_product_mapper[id]   
            #print('after: ', elements)    #debug
            print(*elements, sep='\t')
        else:
            #print('---') #debug
            print(*elements, sep='\t')


