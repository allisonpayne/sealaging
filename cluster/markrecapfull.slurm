#!/bin/bash

#SBATCH --job-name=MRfull
#SBATCH --mail-user=mczapans@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

module load miniconda3.9
conda activate sealaging

datetime=$(date +"%Y%m%d%H%M")
Rscript cluster/markrecapture.R full cluster/output$datetime

echo Transfer outputs locally with:
echo scp -r mczapans@hb.ucsc.edu:sealaging/cluster/output$datetime/ \<path-to-output\>
