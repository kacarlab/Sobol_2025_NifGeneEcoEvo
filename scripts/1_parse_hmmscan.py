#!/usr/bin/env python
# coding: utf-8

# # parse_hmmscan
# 
# Parse hmmscan output to:
# 1. Annotate best KO number per protein in genome dataset

# Usage:
# ./parse_hmmscan.py


import pandas as pd
import glob
import os


def annotate_proteins(files):
    
    # Iterate over hmmscan output files
    for file in files:
    
        base=os.path.basename(file)
        
        fname = os.path.splitext(base)[0]
		
                      
        # Read single hmmscan output
        df = pd.read_csv(file,
                              delim_whitespace = True,
                              skiprows = 3,
                              skipfooter = 10,
                              usecols = [0,2,4],
                              names = ["ko_number", "protein", "e_value"],
                              dtype = {"ko_number": "object",
                                      "protein": "object",
                                      "e_value": "float64"},
                              engine = "python")
        
        # Clean up dataframe
        df["genome"] = df["protein"].str.extract(r"(.+)__")
        df = df[["genome", "protein", "ko_number", "e_value"]]
        
        # Extract rows with best KO number assignment
        # hmmscan output is already sorted by E-value for each protein,
        # so best assignments can be extracted by removing duplicates in protein column
        df.drop_duplicates(subset = ["protein"],
                                keep = "first",
                                inplace = True)

        df["e_value"] = df["e_value"].apply(lambda x: "{:.2e}".format(x))
        
        print(f"Finished parsing {file}")
		
		# Write out protein annotation file
        df.to_csv(f"{fname}.annotated.tsv",
            sep = "\t",
            header = True,
            index = False)	
	
	
# Get hmmscan output file paths
path = "hmmscan_out/"
files = glob.glob(f"{path}*.tbl.out")

print(f"Parsing hmmscan output across {len(files)} files.")

# Annotate
print("Annotating proteins...")

annotate_proteins(files)

print(f"Finished annotating proteins.")    