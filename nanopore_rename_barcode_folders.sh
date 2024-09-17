#!/bin/bash

# Renames barcoded folders by supplying a samplesheet file in the format:
# barcode01     sampleID_xx
# barcode02     sampleID_xx


# Check if the samplesheet argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <samplesheet>"
    exit 1
fi

# Assign the samplesheet argument to a variable
samplesheet="$1"

# Loop through the samplesheet file
while read -r barcode sample_name; do
    # Skip the first two lines containing "library kit" and "flowcell"
    if [[ $barcode == "library" ]] || [[ $barcode == "flowcell" ]]; then
        continue
    fi
    
    # Find the directory with the current barcode name
    for dir in ${barcode}*; do
        # Check if it is a directory
        if [[ -d $dir ]]; then
            # Determine the new directory name
            new_dir="${dir/$barcode/$sample_name}"
            
            # Rename the directory only if the new name is different
            if [[ "$dir" != "$new_dir" ]]; then
                mv "$dir" "$new_dir"
                echo "Renamed directory $dir to $new_dir"
            fi
        fi
    done
done < <(tr -d '\r' < "$samplesheet")
