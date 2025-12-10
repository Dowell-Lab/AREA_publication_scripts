#!/bin/bash 
#SBATCH --job-name=area_run # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ispe1418@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --partition long
#SBATCH --mem=500gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/ispe1418/AREA/eando/area_test.%j.out # Standard output
#SBATCH --error=/scratch/Users/ispe1418/AREA/eando/area_test.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#turn on the virtual machine you are using if you are using one#
path_to_venv=$HOME/VENV/area_run
source $path_to_venv/bin/activate


#set paths to AREA and to files to load in
path_to_area=$HOME/AREA/src/
indir=/Shares/down/public/forAREA/Neuroblastoma/cangelosi_neuroblastoma/
commoncolumn=Participant
rank_file=${indir}NB_rankfile.csv
boolean_attribute_file=${indir}NB_boolean.csv
outdirname=/scratch/Users/ispe1418/AREA/output

mkdir -p "$outdirname"

echo $rank_file
echo $boolean_attribute_file
echo $outdirname

echo $outdirname

python3 ${path_to_area}AREA_core.py --verbose -od $outdirname -cc $commoncolumn -rf $rank_file -baf $boolean_attribute_file --processes 64

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

