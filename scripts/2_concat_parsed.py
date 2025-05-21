import pandas as pd
import os

# Create an empty DataFrame to hold the concatenated data
concatenated_df = pd.DataFrame()

# Directory containing your TSV files
directory = 'hmmscan_parsed/'

# List all files in the directory
file_list = os.listdir(directory)

# Flag to track whether headers have been added
headers_added = False

# Iterate through each file
for file_name in file_list:
    if file_name.endswith('.tsv'):
        file_path = os.path.join(directory, file_name)
        
        # Read the file into a temporary DataFrame
        temp_df = pd.read_csv(file_path, sep='\t')
        
        # Check if headers have been added
        if not headers_added:
            concatenated_df = pd.concat([concatenated_df, temp_df])
            headers_added = True
        else:
            concatenated_df = pd.concat([concatenated_df, temp_df.iloc[1:]])

# Reset the index of the concatenated DataFrame
concatenated_df.reset_index(drop=True, inplace=True)

# Save the concatenated DataFrame to a new TSV file
concatenated_df.to_csv('proteindb.annotated.tsv', sep='\t', index=False)

