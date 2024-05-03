#!/bin/bash

#SBATCH --job-name=MRtest
#SBATCH --mail-user=alrpayne@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --n=4

module load miniconda3.9
conda activate sealaging
Rscript cluster/markrecapture.R full 
