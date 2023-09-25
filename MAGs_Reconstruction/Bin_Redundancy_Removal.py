#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 20:41:38 2023

@author: vishnu
#Same Script was used with replacing 'Cluster 90' with other clustering criteria (eg: 'Cluster 95', 'Cluster 97' and 'Cluster 99') for selecting cluster-representatives
"""
import pandas as pd
import numpy as np
output_file = 'Cluster90_represetatives.txt'
df = pd.read_excel("Final_Bins_Stats_24_09-2023.xlsx", "Sheet1", index_col=0, na_values=["NA"])
# Get the minimum value in 'Column1'
min_value = df['Cluster 90'].min()

# Get the maximum value in 'Column1'
max_value = df['Cluster 90'].max()

with open(output_file, 'a') as file:
    # Get the name of the index column, or provide a custom name if it's not set
    index_column_name = df.index.name if df.index.name is not None else 'Index'

    # Extract the header (column names and index column) as a string
    header_string = '\t'.join([index_column_name] + list(df.columns))
    file.write(header_string + '\n')
    for i in range(1, max_value+1):
        df_2 = df[df["Cluster 90"].isin([i])]
        #print(df_2)
    
        # Find the row with the highest value in Column: completeness
        max_value = df_2['completeness'].max()

        # Filter the DataFrame for rows that meet both criteria
        filtered_df1 = df_2[(df_2['completeness'] == max_value)]
        #print(filtered_df1)

        # Print the number of rows that meet the criteria
        num_rows = len(filtered_df1)
        #print(f"Number of rows that meet the criteria: {num_rows}")
        
        
        # Find the row with the lowest value in Column: contamination
        min_value = filtered_df1['contamination'].min()
        
        if num_rows > 1:
            # Get the row with the maximum value in column 3
            #print(filtered_df1)
            filtered_df2 = filtered_df1[(filtered_df1['contamination'] == min_value)]
            num_rows2 = len(filtered_df2)
            if num_rows2 > 1:
                #print(filtered_df2)
                max_value_N50 = filtered_df2['N50'].max()
                filtered_df3 = filtered_df2[(filtered_df2['N50'] == max_value_N50)]
                print(filtered_df3)
                file.write(filtered_df3.to_csv(sep='\t', index=True, header=False))
            else:
                file.write(filtered_df2.to_csv(sep='\t', index=True, header=False))
        else:
            #print(filtered_df)
            file.write(filtered_df1.to_csv(sep='\t', index=True, header=False))
