#!/bin/bash

# Path to your text file containing tissue names
tpmorreads="reads"
#tpmorreads="tpm"
#$atissue $tpmorreads

indir=/Shares/down/public/omics_databases/GTEX_hg38_RNA/RNA_original_files/
outdir=/Shares/down/public/omics_databases/GTEX_hg38_RNA/RNA_per_tissue/
csv_file="/Shares/down/public/omics_databases/GTEX_hg38_RNA/metadata/tissuesdf_mingreaterthan9samples.csv"

# Skip header and loop through each line
tail -n +2 "$csv_file" | while IFS=',' read -r idx tissue donors prefix; do
    # Remove quotes if present
    prefix_clean=$(echo "$prefix" | tr -d '"')
    
    # Use the variable
    echo "Processing file with prefix: $prefix_clean"
    sbatch --export=prefixfilename=$prefix_clean,tpmorreads=$tpmorreads,outdir=$outdir,indir=$indir outputissue.sbatch
    # Example: do something with $prefix_clean
    # e.g., filename="${prefix_clean}.txt"
done
