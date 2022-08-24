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

# Things to do manually:

# currently nothing

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a "*.gff" file')
#parser.add_argument('-o', '--output', help='path to an output file (default = inputfilename.fixed.gff)')
args = parser.parse_args()

if args.input:
    input_file = Path(args.input)
else:
    print("No input file provided, use '-i' and supply a .gff file")
    input_file = "ArME14_pante.gff3"
    #raise SystemExit

#if args.output:
#    output_file = Path(args.output)
#else:
#    output_file = Path(os.path.join(os.getcwd(),os.path.basename(input_file)+'.fixed.gff'))

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


featureStart = {}
featureStop = {}

previousLine = []
telomere_motives = ['AGGGTT','TAGGGT','TTAGGG','GTTAGG','GGTTAG','GGGTTA','CCCTAA','ACCCTA','AACCCT','TAACCC','CTAACC','CCTAAC']

with open(input_file) as file:
    input = csv.reader(file, delimiter='\t')
    for line in input:
        currentLine = line
        if str(line[0]).startswith('##sequence-region'):
            newLine = str(line[0]).split()
            featureStart[newLine[1]] = newLine[2]
            featureStop[newLine[1]] = newLine[3]
        if currentLine == previousLine:
            continue

        if len(line) > 6:
            # limiting what to look for
            if line[1] == "RepeatMasker" or line[1] == "EAhelitron" or line[1] == "RepeatModeler" or line[1] == "MiteFinderII" or line[1] == "LTRharvest" or line[1] == "LTRdigest" or line[1] == "tRNAScan-SE" or line[1] == "RNAmmer" or line[1] == "pante_protein_families":
                
                elements = []
                for element in line:
                    elements.append(element)

                # MiteFinderII
                if elements[1] == "MiteFinderII":

                    if elements[2] == "terminal_inverted_repeat":
                        elements[2] = "repeat_region"
                        elements[8] = "Name=TIR; rpt_type=inverted; Ontology_term=SO:0000481, SO:terminal_inverted_repeat"
                        #print(*elements, sep='\t')

                    elif elements[2] == "repeat_region":
                        elements[2] = "mobile_element"
                        elements[8] = "Name=MITE; mobile_element_type=transposon:MITE; Ontology_term=SO:0000338, SO:MITE"
                    
                    elif elements[2] == "MITE":
                        continue

                    print(*elements, sep='\t')


                #pante_protein_matches
                if elements[1] == "pante_protein_families":
                    search = re.search('Name=(.+?)$', elements[8])
                    name=search.group(1).split(';')[0]

                    if name == "DNA/TIR/Tc1-Mariner/Fot1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:TcMar-Fot1 "

                    elif name == "LINE/Tad1-like":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:Tad1 "
                        
                    elif name == "LTR/Ty3/Gypsy":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=retrotransposon:Gypsy "

                    elif name == "RC/Helitron":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Helitron-1 "

                    elif name == "DNA/TIR/Tc1-Mariner/Mariner":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:TcMar-Mariner "

                    elif name == "LTR/Ty1/Copia":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=retrotransposon:Copia "

                    elif name == "DNA/TIR/hAT/Hobo":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Hobo "

                    elif name == "DNA/TIR/hAT/Ac":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Activator "

                    elif name == "DNA/TIR/hAT/Tol2":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Tol2 "

                    elif name == "LINE/I":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:L1 "

                    elif name == "LINE/LOA":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:R1-LOA "

                    elif name == "DNA/TIR/MuDR":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:MULE-MuDR "

                    elif name == "DNA/TIR/Tc1-Mariner/Tc1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:TcMar-Tc1 "

                    elif name == "DNA/TIR/Kolobok/IS4EU":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Kolobok-IS4EU "

                    elif name == "DNA/TIR/En-Spm":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:CMC-EnSpm "

                    elif name == "LINE/L1/Tx1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:L1-Tx1 "

                    elif name == "LINE/R1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:R1 "

                    elif name == "LINE/R2":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:R2 "

                    elif name == "LINE/Tad1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:Tad1 "

                    elif name == "LINE/Telomeric":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE "

                    elif name == "LINE/RTE/BovB":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:BovB "
                                                
                    elif name == "LINE/RTE/X":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:X "
                                                                        
                    elif name == "LINE/RTE/RTE":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:RTE "

                    elif name == "LINE/CR1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:CR1 "

                    elif name == "LINE/L2":
                         elements[2] = "mobile_element"
                         elements[8] += "; mobile_element_type=LINE:L2 "

                    elif name == "LINE/Jockey":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:I-Jockey "

                    elif name == "DNA/TIR/Tc1-Mariner":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:TcMar "

                    elif name == "DNA/TIR/Tc1-Mariner/Mogwai":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:TcMar-Mogwai "

                    elif name == "LINE/L1":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:L1 "
                        
                    elif name == "LINE/Proto2":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=LINE:Proto2 "
                                                
                    elif name == "PLE/Penelope":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=retrotransposon:Penelope "

                    elif name == "DNA/TIR/hAT/Restless":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:Restless "

                    elif name == "DNA/TIR/P":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:P-Element "

                    elif name == "DNA/TIR/PIF/Harbinger":
                        elements[2] = "mobile_element"
                        elements[8] += "; mobile_element_type=transposon:PIF-like-Element Harbinger "

                    else:
                                print("#######################################################################")
                                print("Please fix this script and add filters for the before mentioned repeats")
                                print("#######################################################################")
                                print('Add this to script under "pante_protein_families": ',name)
                                raise SystemExit
                    print(*elements, sep='\t')

                # RNAmmer rRNA
                if elements[1] == "RNAmmer":
                    if elements[2] == "rRNA_5S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "Name=5S rRNA; ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "Name=5S rRNA; product=5S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_28S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "Name=28S rRNA; ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "Name=28S rRNA; product=28S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_23S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "Name=23S rRNA; ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "Name=23S rRNA; product=23S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_16S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "Name=16S rRNA; ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "Name=16S rRNA; product=16S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    elif elements[2] == "rRNA_18S":
                        rRNA_parent = re.search('Parent=(.*[0-9]+)$', elements[8])
                        rRNA_p = rRNA_parent.group(1)
                        elements[2] = "gene"
                        elements[8] = "Name=18S rRNA; ID=" + rRNA_p
                        print(*elements, sep='\t')
                        elements[2] = "rRNA"
                        elements[8] = "Name=18S rRNA; product=18S ribosomal RNA; Parent=" + rRNA_p
                        print(*elements, sep='\t')
                    

                # microsatellite
                if elements[2] == "microsatellite":
                    rpt_type = "tandem"
                    elements[2] = "repeat_region"
                    if "repeat_unit" in elements[8]:

                        search = re.search('repeat_unit=(.+?)$', elements[8])
                        
                        # removing 'repeat_unit' tag
                        element_8 = elements[8].split(";")
                        idx = 0
                        for i in element_8:
                            idx +=1
                            if i.startswith("repeat_unit"):
                                break
                        del element_8[idx -1]                        
                        element_8_new = ';'.join(map(str,element_8))                        
                        elements[8] = element_8_new

                        # adding 'rpt_type' and 'satellite' tags
                        if search:
                            rpt_unit = search.group(1)
                            # telomere detection
                            if rpt_unit in telomere_motives and int(featureStart[elements[0]] + 1000) > int(elements[3]):
                                elements[8] = "Name=telomere; Ontology_term=SO:0000624, SO:telomere; rpt_type=telomeric_repeat; satellite=microsatellite"
                                # print(elements)

                            elif rpt_unit in telomere_motives and int(featureStop[elements[0]] - 1000 ) < int(elements[4]):
                                elements[8] = "Name=telomere; Ontology_term=SO:0000624, SO:telomere; rpt_type=telomeric_repeat; satellite=microsatellite"
                                # print(elements)


                            elif int(elements[4])-int(elements[3]) < len(rpt_unit):
                                continue
                            else:
                                elements[8] += "; rpt_type=" + rpt_type + "; satellite=microsatellite"
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite=microsatellite"
                    
                    print(*elements, sep='\t')

                
                # minisatellite
                if elements[2] == "minisatellite":
                    rpt_type = "tandem"
                    elements[2] = "repeat_region"
                    if "repeat_unit" in elements[8]:
                        search = re.search('repeat_unit=(.+?)$', elements[8])

                        # removing 'repeat_unit' tag
                        element_8 = elements[8].split(";")
                        idx = 0
                        for i in element_8:
                            idx +=1
                            if i.startswith("repeat_unit"):
                                break
                        del element_8[idx -1]                        
                        element_8_new = ';'.join(map(str,element_8))                        
                        elements[8] = element_8_new
                        # adding 'rpt_type' and 'satellite' tags

                        if search:
                            rpt_unit = search.group(1)
                            elements[8] += "; rpt_type=" + rpt_type + "; satellite=minisatellite"
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite=minisatellite"
                    
                    print(*elements, sep='\t')


                # monomeric_repeat
                if elements[2] == "monomeric_repeat":
                    rpt_type = "tandem"
                    elements[2] = "repeat_region"
                    if "repeat_unit" in elements[8]:
                        search = re.search('repeat_unit=(.+?)$', elements[8])

                        # removing 'repeat_unit' tag
                        element_8 = elements[8].split(";")
                        idx = 0
                        for i in element_8:
                            idx +=1
                            if i.startswith("repeat_unit"):
                                break
                        del element_8[idx -1]                        
                        element_8_new = ';'.join(map(str,element_8))                        
                        elements[8] = element_8_new
                        # adding 'rpt_type' and 'satellite' tags
                        if search:
                            rpt_unit = search.group(1)
                            elements[8] += "; rpt_type=" + rpt_type + "; satellite=satellite"
                    else:
                        elements[8] += "; rpt_type=" + rpt_type + "; satellite=satellite"

                    print(*elements, sep='\t')

                # tRNA
                if elements[1] == "tRNAScan-SE":
                                        
                    trna = re.search('(.*)_tRNA', elements[2]) 
                    if trna:
                        x = trna.group()
                        if x in trnas:
                            parent = re.search('Parent=.*[0-9]+', elements[8])
                            p = parent.group(0)

                            #to overwrite the anticodon tag
                            #anticodon_tag = re.search('anticodon=[A-Z]{3};', elements[8])
                            #anticodon_tag_old = anticodon_tag.group()
                            #elements[8] = elements[8].replace(anticodon_tag_old, "anticodon" + anticodon_tag)
                            
                            if 'pseudo' in p:
                                elements[2] = "gene"
                                gene = trnas[x][1]
                                id = p.split("=",1)[1]
                                id1 = id.split("_",1)[0]
                                elements[8] = "ID=trna_gene." + id1 + "; gene=" + gene + "; pseudo=true; pseudogene=unknown"
                                print(*elements, sep='\t')
                                elements[2] = "trna"
                                product = trnas[x][0]
                                elements[8] = "ID=trna." + id1 + "; Parent=gene." + id1 + ";product=" + product + "; pseudo=true; pseudogene=unknown"
                                print(*elements, sep='\t')
                                
                            else:
                                elements[2] = "gene"
                                gene = trnas[x][1]
                                id = p.split("=",1)[1]
                                id1 = id.split("_",1)[1]
                                elements[8] = "ID=trna_gene." + id1 + "; gene=" + gene
                                print(*elements, sep='\t')
                                elements[2] = "trna"
                                product = trnas[x][0]
                                elements[8] = "ID=trna." + id1 + "; Parent=gene." + id1 + ";product=" + product
                                print(*elements, sep='\t')
                                

                # helitron
                if elements[2] == "helitron":
                    elements[2] = "mobile_element"
                    # removing 'ID' tag
                    element_8 = elements[8].split(";")
                    idx = 0
                    for i in element_8:
                        idx +=1
                        if i.startswith("ID="):
                            break
                    del element_8[idx -1]                        
                    element_8_new = ';'.join(map(str,element_8))                        
                    elements[8] = element_8_new
                    element_8 = elements[8].split(";")
                    idx = 0
                    for i in element_8:
                        idx +=1
                        if i.startswith("Parent="):
                            break
                    del element_8[idx -1]                        
                    element_8_new = ';'.join(map(str,element_8))                        
                    elements[8] = element_8_new

                    elements[8] = "Name=helitron; Ontology_term=SO:0000544, SO:helitron; mobile_element_type=transposon:Helitron"

                    print(*elements, sep='\t')

               # LTR harvest
               
                # LTRharvest
                if elements[2] == "LTR_retrotransposon":
                    elements[2] = "mobile_element"
                    elements[8] = "Name=LTR retrotransposon; mobile_element_type=retrotransposon; Ontology_term=SO:0000186"
                    print(*elements, sep='\t')
                
                if elements[2] == "long_terminal_repeat":
                    elements[2] = "repeat_region"
                    elements[8] = "Name=LTR; rpt_type=long_terminal_repeat; rpt_family=LTR; Ontology_term=SO:0000286"
                    print(*elements, sep='\t')

                if elements[2] == "target_site_duplication":
                    elements[2] = "repeat_region"
                    elements[8] = "Name=Target site duplication; rpt_type=direct; rpt_family=LTR; Ontology_term=SO:0000434"
                    print(*elements, sep='\t')

                # repeat_region
                if elements[2] == "repeat_region":
                    if "repeat_family" in elements[8]:
                        search = re.search('repeat_family=(.+?)$', elements[8])
                        if search:
                            rpt_family = search.group().split(";")[0].split("=")[1]
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
                                elements[8] += "; mobile_element_type=transposon:Activator "

                            elif rpt_family.lower() == "dna/hat-charlie":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Charlie "

                            elif rpt_family.lower() == "dna/hat-tip100":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Tip100 "

                            elif rpt_family.lower() == "dna/ginger":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Ginger "

                            elif rpt_family.lower() == "dna/ginger-2":
                                elements[2] = "mobile_element"
                                elements[8] += "; mobile_element_type=transposon:Ginger-2 "

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
                                elements[8] += "; rpt_type=other; note=potential rRNA  "

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
                                print('Add this to script under "repeat_region": ',rpt_family)
                                raise SystemExit

                        print(*elements, sep='\t')
        previousLine = currentLine
            #else:
                #print(line[1])


                
                # print the modified list of repeat elements with a tab as the separator
                #print(*elements, sep='\t')
    

