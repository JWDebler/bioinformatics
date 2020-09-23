import os
import csv
import argparse
from pathlib import Path

# This script extracts repeat annotations made via the PiRATE pipeline from a gff file
# and creates a gff file that can be used by table2asn for annotation submissions
# to genbank. Definitions are taken from https://www.ncbi.nlm.nih.gov/WebSub/html/annot_examples.html
# but this does not seem like an exhaustive list.

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
parser.add_argument('-o', '--output', help='path to an output file (default = inputfilename.fixed.gff)')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    print("No input file provided, use '-i' and supply a .gff file")
    raise SystemExit

if args.output:
    output_file = Path(args.output)
else:
    output_file = Path(os.path.join(os.getcwd(),os.path.basename(input_file)+'.fixed.gff'))

gff = {}


#If file with same outputname exists, delete it first
if os.path.isfile(output_file):
    os.unlink(output_file)


with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        #look for the 'PiRATE' identifier and then check for the type of repeat
        if line[1] == "PiRATE":
            elements = []
            for element in line:
                elements.append(element)

            if elements[2] == "LTR":
                rpt_type = "long_terminal_repeat"
                elements[2] = "repeat_region"
                if "Gypsy" in line[8]:
                    rpt_family = "Gypsy"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family
                    
                elif "Copia" in line[8]:
                    rpt_family = "Copia"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family
                    
                else:
                    elements[8]+="; rpt_type="+rpt_type
                                        

            elif elements[2] == "noCat":
                rpt_type = "other"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type="+rpt_type

            elif elements[2] == "LTR|TIR":
                rpt_type = "other"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type="+rpt_type

            elif elements[2] == "TRIM":
                rpt_type = "long_terminal_repeat"
                rpt_family = "TRIM"
                elements[2] = "repeat_region"
                elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

            elif elements[2] == "DIRS":
                rpt_type = "long_terminal_repeat"
                rpt_family = "DIRS"
                elements[2] = "repeat_region"
                elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

            elif elements[2] == "LINE":
                rpt_type = "other"
                rpt_family = "LINE"
                elements[2] = "repeat_region"
                elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

            elif elements[2] == "SINE":
                rpt_type = "dispersed"
                rpt_family = "SINE"
                elements[2] = "repeat_region"
                elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

            elif elements[2] == "MITE":
                rpt_type = "other"
                rpt_family = "transposon: MITE"
                elements[2] = "repeat_region"
                elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

            elif elements[2] == "TIR":
                rpt_type = "terminal"
                elements[2] = "repeat_region"
                if "Tc1-Mariner" in line[8]:
                    rpt_family = "Tc1-Mariner"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

                elif "CACTA" in line[8]:
                    rpt_family = "CACTA"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family
            
                elif "MuDR" in line[8]:
                    rpt_family = "MuDR"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

                elif "PIF-Harbinger" in line[8]:
                    rpt_family = "PIF-Harbinger"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family

                else:
                    elements[8]+="; rpt_type="+rpt_type
                    
            elif elements[2] == "SSR":
                rpt_type = "other"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type="+rpt_type

            elif elements[2] == "LARD":
                rpt_type = "long_terminal_repeat"
                elements[2] = "repeat_region"
                mobile_element = "retrotransposon: LARD"
                elements[8]+="; rpt_type=" + rpt_type + "; mobile_element=" + mobile_element

            elif elements[2] == "Helitron":
                rpt_type = "dispersed"
                mobile_element = "transposon: Helitron"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type=" + rpt_type + "; mobile_element=" + mobile_element

            elif elements[2] == "Maverick":
                rpt_type = "dispersed"
                mobile_element = "transposon: Maverick"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type=" + rpt_type + "; mobile_element=" + mobile_element
            # print the modified list of repeat elements with a tab as the separator
            print(*elements, sep='\t')
        
        else:
            # print all the other elements with a tab as the separator
            print(*line, sep='\t')  

