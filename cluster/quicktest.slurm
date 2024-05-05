#!/bin/bash

#SBATCH --job-name=quicktest
#SBATCH --mail-user=mczapans@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=5
#SBATCH --ntasks=1
#SBATCH --mem=20M
#SBATCH --cpus-per-task=1

module load miniconda3.9
conda activate sealaging

# Output folder
datetime=$(date +"%Y%m%d%H%M")
#Rscript cluster/markrecapture.R test cluster/output$datetime
echo Copy results locally with:
echo scp mczapans@hb.ucsc.edu:sealaging/cluster/output$datetime \<path-to-output/>
