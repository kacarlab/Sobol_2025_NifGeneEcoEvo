#!/bin/bash

# Usage example:
# ./run_hmmer.sh

# Set up conda environment

cp /staging/msobol/conda_pkg/hmmer.tar.gz ./

# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=hmmer
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR # Extract conda environment tar package
. $ENVDIR/bin/activate


# Extract profile database
# Should extract 5 files: kofam, kofam.h3f, kofam.h3i, kofam.h3m, kofam.h3p
cp /staging/msobol/databases/kofam.tar.gz ./ # Copy from staging directory
tar -xzf kofam.tar.gz


# Renaming first argument
file=$1
bname=$(basename "$file" .targetseqs.fasta)

# Run hmmscan

echo "Conducting HMM scan."

hmmscan --tblout ${bname}.tbl.out -E 0.00001 --domE 0.00001 --incE 0.00001 --incdomE 0.00001 --noali kofam ${file}

echo "Finished conducting HMM scan."

cp ${bname}.tbl.out /staging/msobol/outputs/hmmer/

# Clean up
echo "Clean up"

rm -r kofam* hmmer* *fasta *.out *.sh

#END
