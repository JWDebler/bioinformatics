import os
import csv
import argparse
import re
from pathlib import Path

# This script extracts repeat annotations made via the panTE pipeline from a gff file
# and creates a gff file that can be used by table2asn for annotation submissions
# to genbank. Definitions are taken from https://www.ncbi.nlm.nih.gov/WebSub/html/annot_examples.html
# but this does not seem like an exhaustive list.
# http://www.insdc.org/documents/feature-table
# https://github.com/Dfam-consortium/RepeatModeler/blob/master/RepeatClassifier
# http://www.insdc.org/controlled-vocabulary-rpttype-qualifier


#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
parser.add_argument('-o', '--output', help='path to an output file (default = inputfilename.fixed.gff)')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    input_file = ("GCA_004011705.1_Alentis_Al4_genomic_clean_header_pante.gff3")
    print("No input file provided, use '-i' and supply a .gff file")
    #raise SystemExit

if args.output:
    output_file = Path(args.output)
else:
    output_file = Path(os.path.join(os.getcwd(),os.path.basename(input_file)+'.fixed.gff'))

gff = {}
trnas = {"tRNA":["tRNA-Xxx", "trnX"], 
        "alanyl_tRNA": ["tRNA-Ala", "trnA"],
        "glutaminyl_tRNA": ["tRNA-Gln", "trnQ"],
        "prolyl_tRNA": ["tRNA-Pro", "trnP"],
        "glytamyl_tRNA":["tRNA-Glu", "trnE"],
        "methionyl_tRNA": ["tRNA-Met", "trnM"],
        "asparaginyl_tRNA": ["tRNA-Asn", "trnN"],
        "threonyl_tRNA": ["tRNA-Thr", "trnT"],
        "glycyl_tRNA": ["tRNA-Gly","trnG"],
        "valyl_tRNA": ["tRNA-Val", "trnV"],
        "tyrosyl_tRNA": ["tRNA-Tyr","trnY"],
        "cysteinyl_tRNA": ["tRNA-Cys","trnC"],
        "isoleucyl_tRNA": ["tRNA-Ile", "trnI"],
        "seryl_tRNA": ["tRNA-Ser", "trnS"],
        "leucyl_tRNA": ["tRNA-Leu","trnL"],
        "selenocysteinyl_tRNA": ["tRNA-Sec","trnU"],
        "tryptophanyl_tRNA": ["tRNA-Trp","trnW"],
        "pyrrolysyl_tRNA": ["tRNA-Pyl","trnO"],
        "lysyl_tRNA": ["tRNA-Lys","trnK"],
        "aspartyl_tRNA": ["tRNA-Asp","trnD"],
        "arginyl_tRNA": ["tRNA-Arg","trnR"],
        "histidyl_tRNA": ["tRNA-His","trnH"],
        "phenylalanyl_tRNA": ["tRNA-Phe", "trnF"]}

#If file with same outputname exists, delete it first
if os.path.isfile(output_file):
    os.unlink(output_file)


with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
                    

        if len(line) > 6:
            if line[1] == "RepeatMasker" or line[1] == "EAhelitron" or line[1] == "RepeatModeler" or line[1] == "MiteFinderII" or line[1] == "LTRharvest" or line[1] == "LTRdigest" or line[1] == "tRNAScan-SE" or line[1] == "RNAmmer":
                
                elements = []
                for element in line:
                    elements.append(element)

                # RNAmmer rRNA
                if elements[1] == "RNAmmer":
                    if elements[2] == "rRNA_5S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "product=5S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_28S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "product=28S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_23S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "product=23S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_16S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "product=16S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_18S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "product=18S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    

                # microsatellite
                if elements[2] == "microsatellite":
                    rpt_type = "tandem"
                    elements[2] = "microsatellite"
                    if "repeat_unit" in elements[8]:
                        search = re.search('repeat_unit=(.+?)$', elements[8])
                        if search:
                            rpt_unit = search.group(1)
                            if elements[4]-elements[3] < rpt_unit.length():
                                print("ERROR")
                                break
                            elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=microsatellite" + "; rpt_unit_seq=" + rpt_unit
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=microsatellite"
                    
                    print(*elements, sep='\t')

                # minisatellite
                if elements[2] == "minisatellite":
                    rpt_type = "tandem"
                    elements[2] = "minisatellite"
                    if "repeat_unit" in elements[8]:
                        search = re.search('repeat_unit=(.+?)$', elements[8])
                        if search:
                            rpt_unit = search.group(1)
                            elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=minisatellite" + "; rpt_unit_seq=" + rpt_unit
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=minisatellite"
                    
                    print(*elements, sep='\t')


                # monomeric_repeat
                if elements[2] == "monomeric_repeat":
                    rpt_type = "tandem"
                    elements[2] = "repeat_region"
                    if "repeat_unit" in elements[8]:
                        search = re.search('repeat_unit=(.+?)$', elements[8])
                        if search:
                            rpt_unit = search.group(1)
                            elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=satellite" + "; rpt_unit_seq=" + rpt_unit
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite_type=satellite"

                    print(*elements, sep='\t')

                # tRNA
                if elements[1] == "tRNAScan-SE":
                                        
                    trna = re.search('(.*)_tRNA', elements[2]) 
                    if trna:
                        x = trna.group()
                        if x in trnas:
                            parent = re.search('Parent=.*[0-9]+', elements[8])
                            p = parent.group(0)
                            anticodon_tag = re.search('anticodon=[A-Z]{3};', elements[8])
                            anticodon_tag_old = anticodon_tag.group()
                            elements[8] = elements[8].replace(anticodon_tag_old, "")
                            if 'pseudo' in p:
                                elements[2] = "trna"
                                product = trnas[x][0]
                                elements[8] += "; product=" + product + "; pseudo=true; pseudogene=unknown"
                                print(*elements, sep='\t')
                                elements[2] = "gene"
                                gene = trnas[x][1]
                                id = p.split("=",1)[1]
                                elements[8] = "ID=" + id + "; gene=" + gene + "; pseudo=true; pseudogene=unknown"
                                print(*elements, sep='\t')
                            else:
                                elements[2] = "trna"
                                product = trnas[x][0]
                                elements[8] += "; product=" + product
                                print(*elements, sep='\t')
                                elements[2] = "gene"
                                gene = trnas[x][1]
                                id = p.split("=",1)[1]
                                elements[8] = "ID=" + id + "; gene=" + gene
                                print(*elements, sep='\t')

                # helitron
                if elements[2] == "helitron":
                    elements[2] = "mobile_element"
                    elements[8] += "; mobile_element_type=transposon:Helitron"

                    print(*elements, sep='\t')

               # LTR harvest
                if elements[2] == "LTR_retrotransposon":
                    elements[2] = "mobile_element"
                    elements[8] = "rpt_type=long_terminal_repeat; mobile_element_type=retrotransposon; rpt_family=LTR"

                    print(*elements, sep='\t')

                # repeat_region
                if elements[2] == "repeat_region":
                    if "repeat_family" in elements[8]:
                        search = re.search('repeat_family=(.+?)$', elements[8])
                        if search:
                            rpt_family = search.group(1)
                            if rpt_family.lower() == "ltr/gypsy":
                                elements[2] =  "mobile_element"
                                elements[8] += "; mobile_element_type=retrotransposon:Gypsy"

                            elif rpt_family.lower() == "ltr/copia":
                                elements[2] =  "mobile_element"
                                elements[8] += "; mobile_element_type=retrotransposon:Copia"

                            elif rpt_family.lower() == "ltr":
                                elements[8] += "; rpt_type=long_terminal_repeat; rpt_family=LTR"

                            elif rpt_family.lower() == "dna/tcmar-fot1":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:TcMar-Fot1"

                            elif rpt_family.lower() == "dna/tcmar-tc1":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:TcMar-Tc1"

                            elif rpt_family.lower() == "dna/tcmar-tc2":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:TcMar-Tc2"

                            elif rpt_family.lower() == "dna/tcmar-tc4":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:TcMar-Tc4"

                            elif rpt_family.lower() == "dna/cmc-enspm":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:CMC-EnSpm "

                            elif rpt_family.lower() == "dna/mule-mudr":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:MULE-MuDR "

                            elif rpt_family.lower() == "dna/merlin":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Merlin  "

                            elif rpt_family.lower() == "dna/hat-ac":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:hAT-Ac "

                            elif rpt_family.lower() == "line/penelope":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=LINE:penelope "

                            elif rpt_family.lower() == "line/i-jockey":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=LINE:I-Jockey "

                            elif rpt_family.lower() == "line/l1":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=LINE:L1 "

                            elif rpt_family.lower() == "line/l2":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=LINE:L2 "

                            elif rpt_family.lower() == "line/tad1":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=LINE:Tad1 "

                            elif rpt_family.lower() == "artefact":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=other:Artifact  "

                            elif rpt_family.lower() == "rrna":
                                elements[2] = "repeat_region"
                                elements[8] += "; rpt_type=other  "

                            elif rpt_family.lower() == "unknown":
                                elements[2] = "repeat_region"
                                elements[8] += "; rpt_type=dispersed  "

                            elif rpt_family.lower() == "dna":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:DNA  "

                            elif rpt_family.lower() == "dna/kolobok":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Kolobok  "

                            elif rpt_family.lower() == "dna/kolobok-e":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Kolobok-E  "

                            elif rpt_family.lower() == "dna/kolobok-h":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Kolobok-H  "

                            elif rpt_family.lower() == "dna/kolobok-hydra":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Kolobok-Hydra-specific_Branch  "

                            elif rpt_family.lower() == "dna/kolobok-t2":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Kolobok-T2 "
                            
                            elif rpt_family.lower() == "rc/helitron":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Helitron-1 "

                            elif rpt_family.lower() == "rc/helitron-2":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Helitron-2 "

                            else:
                                print("#######################################################################")
                                print("Please fix this script and add filters for the before mentioned repeats")
                                print("#######################################################################")
                                print('Add this to script: ',rpt_family)
                                raise SystemExit

                        print(*elements, sep='\t')

                        print("deeeeeeeeeeerp")
                        
            #else:
                #print(line[1])


                
                # print the modified list of repeat elements with a tab as the separator
                #print(*elements, sep='\t')
    

