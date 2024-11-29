#!/bin/bash

# Check if a log file is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <log_file>"
  exit 1
fi

# Use the first argument as the log file
log_file="$1"

# Define start of table
start_line=$(grep -n "Summary of Genetic Correlation Results" "$log_file" | cut -d: -f1)

# Extract all SNP counts
snps_counts=($(grep -o "[0-9]\+ SNPs with valid alleles\." "$log_file" | grep -o "[0-9]\+"))

# Extract table, add SNP count as last column, and print to stdout
awk -v start="$start_line" -v snps="${snps_counts[*]}" '
    BEGIN {
        # Split the SNP counts string into an array
        split(snps, snp_array, " ")
        snp_idx=1
    }
    NR > start {
        # Stop if we encounter a blank line or "Analysis finished at"
        if (/^[[:space:]]*$/ || /Analysis finished at/) exit

        # Add "NrSnps" header to the first line only
        if (NR == start + 1) {
            print $0 " NrSnps"
        } else if (NF) {
            # Print lines in the table with the SNP count appended
            print $0 " " snp_array[snp_idx];
            snp_idx++
        }

    }
' "$log_file" 
