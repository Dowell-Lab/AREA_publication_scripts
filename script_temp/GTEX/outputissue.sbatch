#!/bin/bash
#SBATCH --job-name=tissuesummary
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=05:00:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=50gb                                       # Memory limit
#SBATCH --output=/scratch/Users/allenma/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/allenma/e_and_o/%x_%j.err

################################################################                                                                                                                      
##Load required modules 
################################################################
module load R/4.4.0


################################################################ 
##Set paths
################################################################
rscript_path=/Shares/down/public/omics_databases/GTEX_hg38_RNA/scripts/outputissue.R
################################################################ 
##Job information
################################################################ 
printf "\nR code being run: ${rscript_path}"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou have requested $SLURM_CPUS_ON_NODE core(s).\n"
printf "\n$atissue\n"
printf "\n$tpmorreads\n"

################################################################                  
## Run Script 
################################################################
#using the Rscript command to run the script
Rscript ${rscript_path} "$prefixfilename" $tpmorreads $indir $outdir

echo Job finished at `date +"%T %a %d %b %Y"`


