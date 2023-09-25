#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 13:59:09 2023

@author: vishnu
"""
import pandas as pd
import numpy as np
import os
import glob

# Define the directory path
directory_path = '/Users/vishnu/WORK/METAGENOMIC_DATA_ANALYSIS/FIBER_METAGENOMIC_DATA_ANALYSIS/08_SAMPLEWISE_BINNING/FINAL_BINS'

# Use glob to list only .fasta files in the directory
file_list = glob.glob(os.path.join(directory_path, '*_represetatives.txt'))

# Loop through the .fasta files and open each one
for file_path in file_list:
    #output_file = 'Cluster_Summary.txt'
    output_file = f"{file_path}_summary.tsv"
    # Open the file
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    hi_complete = len(df[(df['completeness'] >= 95) & (df['contamination'] <= 5)])

    ne_complete = len(df[(df['completeness'] >= 90) & (df['completeness'] < 95) & (df['contamination'] <= 10)])
    me_complete = len(df[(df['completeness'] >= 70) & (df['completeness'] < 90) & (df['contamination'] <= 10)])
    partial = len(df[(df['completeness'] < 70) & (df['contamination'] <= 10)])
    
    # Open the output file in write mode ('w') using a 'with' statement
    with open(output_file, 'w') as file:
        # Write the content to the file
        file.write(f"Number of highly complete bins: {hi_complete}\n")
        file.write(f"Number of nearly complete bins: {ne_complete}\n")
        file.write(f"Number of medium complete bins: {me_complete}\n")
        file.write(f"Number of partial bins: {partial}\n")
