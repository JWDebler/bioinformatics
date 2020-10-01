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
    input_file = ("Alentis_AllFeatures_GFF_correct_repeats.gff")
    #print("No input file provided, use '-i' and supply a .gff file")
    #raise SystemExit

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
                if "Gypsy" in line[8]:
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=retrotransposon:Gypsy"
                elif "Copia" in line[8]:
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=retrotransposon:Copia"
                else:
                    elements[2] = "repeat_region"
                    elements[8] += "; rpt_type=long_terminal_repeat; rpt_family=LTR"
                                        

            elif elements[2] == "noCat":
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type=other"

            elif elements[2] == "LTR|TIR":
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type=other"

            elif elements[2] == "TRIM":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=retrotransposon:TRIM"
                
            elif elements[2] == "DIRS":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=retrotransposon:DIRS"
                
            elif elements[2] == "LINE":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=LINE "

            elif elements[2] == "SINE":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=SINE "

            elif elements[2] == "MITE":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=MITE "

            elif elements[2] == "TIR":
                rpt_type = "terminal"
                elements[2] = "repeat_region"

                if "Tc1-Mariner" in line[8]:
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=transposon:TcMar-Tc1"

                elif "CACTA" in line[8]:
                    rpt_family = "CACTA"
                    elements[8] += "; rpt_type=" + rpt_type + "; rpt_family=" + rpt_family
            
                elif "MuDR" in line[8]:
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=transposon:MULE-MuDR "

                elif "PIF-Harbinger" in line[8]:
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=transposon:PIF_Harbinger "

                else:
                    elements[8]+="; rpt_type="+rpt_type
                    
            elif elements[2] == "SSR":
                rpt_type = "other"
                elements[2] = "repeat_region"
                elements[8]+="; rpt_type="+rpt_type

            elif elements[2] == "LARD":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=retrotransposon:LARD"
            
            elif elements[2] == "Helitron":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=transposon:Helitron"

            elif elements[2] == "Maverick":
                elements[2] = "mobile_element"
                elements[8] += "; mobile_element_type=transposon:Maverick"
            
            else:
                print("########################################################################")
                print("# Please fix this script and add filters for the following repeat type #")
                print("########################################################################")
                print(line)
                print("Look for this: " + line[2])
                print("########################################################################")
                raise SystemExit

            print(*elements, sep='\t')
        
        #else:
            # print all the other elements with a tab as the separator
            # print(*line, sep='\t')  

