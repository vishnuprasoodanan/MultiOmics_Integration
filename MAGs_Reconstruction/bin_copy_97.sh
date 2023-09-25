#!/bin/bash

# Source folder
source_folder="/Users/vishnu/WORK/METAGENOMIC_DATA_ANALYSIS/FIBER_METAGENOMIC_DATA_ANALYSIS/08_SAMPLEWISE_BINNING/CLUSTERING_BINS/ALL_BINS"

# Destination folder
destination_folder="/Users/vishnu/WORK/METAGENOMIC_DATA_ANALYSIS/FIBER_METAGENOMIC_DATA_ANALYSIS/08_SAMPLEWISE_BINNING/CLUSTERING_BINS/ALL_BINS/CLUSTER_97_Representatives"

# Text file containing file names to copy (one file per line)
file_list="/Users/vishnu/WORK/METAGENOMIC_DATA_ANALYSIS/FIBER_METAGENOMIC_DATA_ANALYSIS/08_SAMPLEWISE_BINNING/CLUSTERING_BINS/ALL_BINS/Cluster97_represetatives.txt_list.txt"

# Loop through each line in the text file
while IFS= read -r file_name; do
    # Check if the file exists in the source folder
    if [ -e "$source_folder/$file_name" ]; then
        # Copy the file to the destination folder
        cp "$source_folder/$file_name" "$destination_folder/"
        echo "Copied: $file_name"
    else
        echo "File not found: $file_name"
    fi
done < "$file_list"

