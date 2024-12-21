#!/bin/bash

# Renames barcoded files by supplying a samplesheet file in the format:
# barcode01     sampleID_xx
# barcode02     sampleID_xx

# Check if the samplesheet file is provided as an argument
if [[ -z "$1" ]]; then
    echo "Usage: $0 <samplesheet>"
    exit 1
fi

# Define the path to your samplesheet
samplesheet="$1"

# Loop through the samplesheet file
while read -r barcode sample_name; do
    # Skip lines that do not start with 'barcode'
    if [[ ! $barcode =~ ^barcode ]]; then
        continue
    fi
    
    # Find all files starting with the current barcode
    for file in ${barcode}.*; do
        # Ensure $file is a regular file and not a directory
        if [[ -f $file ]]; then
            # Rename the file by replacing the barcode with the sample name
            mv "$file" "${file/$barcode/$sample_name}"
            echo "Renamed $file to ${file/$barcode/$sample_name}"
        fi
    done
done < "$samplesheet"
