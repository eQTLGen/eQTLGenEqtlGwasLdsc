#!/bin/bash

# Check if a log file is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <log_file>"
  exit 1
fi


# Use the first argument as the log file
log_file="$1"

# Extract the block of text for "Heritability of phenotype 1" only
phenotype_block=$(awk '/Heritability of phenotype 1/,/Heritability of phenotype 2\/[0-9]+/{
    if ($0 ~ /Heritability of phenotype 2\/[0-9]+/) exit
    print
}' "$log_file")

# Extract each field and assign to variables
h2=$(echo "$phenotype_block" | grep "Total Observed scale h2" | awk '{print $5}')
h2_se=$(echo "$phenotype_block" | grep "Total Observed scale h2" | awk -F'[()]' '{print $2}')
lambda=$(echo "$phenotype_block" | grep "Lambda GC" | awk '{print $3}')
meanChi2=$(echo "$phenotype_block" | grep "Mean Chi\^2" | awk '{print $3}')
intercept=$(echo "$phenotype_block" | grep "Intercept:" | awk '{print $2}')
intercept_se=$(echo "$phenotype_block" | grep "Intercept:" | awk -F'[()]' '{print $2}')
ratio=$(echo "$phenotype_block" | grep "Ratio" | awk '{print $2}')
ratio_se=$(echo "$phenotype_block" | grep "Ratio" | awk -F'[()]' '{print $2}')

# Print header
echo "log_file	h2	h2_se	lambda	meanChi2	intercept	intercept_se	ratio	ratio_se"

# Print extracted values in tab-separated format
echo "$log_file	$h2	$h2_se	$lambda	$meanChi2	$intercept	$intercept_se	$ratio	$ratio_se"
