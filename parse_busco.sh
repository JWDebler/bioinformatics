#!/bin/bash
# Pareses the BUSCO output file and prints results to stout

# Check if input file is supplied as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Input file from the command line argument
input_file="$1"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File $input_file not found!"
    exit 1
fi

# Extract the sample name (part before the first ".")
sample_name=$(basename "$input_file" | cut -d'.' -f1)

# Extract the relevant line that contains the percentage values
result_line=$(grep -oP 'C:\K[^\n]*' "$input_file")

# Use regular expressions to extract values for each category
complete=$(echo "$result_line" | grep -oP '[0-9.]+(?=%\[S)')
complete_single=$(echo "$result_line" | grep -oP '(?<=S:)[0-9.]+(?=%)')
complete_duplicated=$(echo "$result_line" | grep -oP '(?<=D:)[0-9.]+(?=%)')
fragmented=$(echo "$result_line" | grep -oP '(?<=F:)[0-9.]+(?=%)')
missing=$(echo "$result_line" | grep -oP '(?<=M:)[0-9.]+(?=%)')

# Output the values to stdout in TSV format
#echo -e "sample\tcomplete\tcomplete single\tcomplete duplicated\tfragmented\tmissing"
echo -e "$sample_name\t$complete\t$complete_single\t$complete_duplicated\t$fragmented\t$missing"