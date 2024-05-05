#!/bin/bash

#SBATCH --job-name=MRtest
#SBATCH --mail-user=mczapans@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=30
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

module load miniconda3.9
conda activate sealaging

datetime=$(date +"%Y%m%d%H%M")
Rscript cluster/markrecapture.R test cluster/output$datetime

echo Transfer outputs locally with:
echo scp mczapans@hb.ucsc.edu:sealaging/
