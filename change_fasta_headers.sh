#!/bin/bash

# Check for arguments
if [ $# -lt 2 ]; then
  echo "Usage: $0 <fasta_file> --species <species_name>"
  exit 1
fi

# Assign variables
fasta_file="$1"
species=""
id=$(basename "$fasta_file" .fasta)  # Extract the filename before .fasta

# Parse the arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --species)
      species="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done

# Check if species is provided
if [ -z "$species" ]; then
  echo "Error: --species argument is required."
  exit 1
fi

# Use a temporary file for in-place modification
tmp_file=$(mktemp)

# Process the FASTA file
awk -v species="$species" -v id="$id" '
  /^>/ {
    if ($0 ~ /^>ctg[0-9]+/) {
      sub(/^>ctg[0-9]+/, "& " species " isolate " id, $0)
    } 
    else if ($0 ~ /^>mito/) {
      sub(/^>mito/, "& " species " isolate " id, $0)
    } 
    else if ($0 ~ /^>plasmid/) {
      sub(/^>plasmid/, "& " species " isolate " id, $0)
    }
    print
    next
  }
  { print }  # Print sequence lines unchanged
' "$fasta_file" > "$tmp_file"

# Overwrite the input file with the modified content
mv "$tmp_file" "$fasta_file"

echo "Headers have been modified in-place: $fasta_file"
