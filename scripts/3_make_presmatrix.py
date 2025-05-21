#!/usr/bin/env python
# coding: utf-8

import pandas as pd

# Read the annotation file
annot_df = pd.read_csv("proteindb.annotated.tsv", sep="\t", dtype=str)

# Ensure the file has the required columns
if "genome" not in annot_df.columns or "ko_number" not in annot_df.columns:
    raise ValueError("The input file must contain 'genome' and 'ko_number' columns.")

# Create the presence/absence matrix
pres_abs_df = annot_df.groupby(["genome", "ko_number"]).size().unstack(fill_value=0)

# Convert counts to binary presence/absence (1 if present, 0 if absent)
pres_abs_df = (pres_abs_df > 0).astype(int)

# Save the presence/absence matrix
pres_abs_df.to_csv("proteindb.presmatrix.tsv", sep="\t")

print("Finished generating presence/absence matrix.")
