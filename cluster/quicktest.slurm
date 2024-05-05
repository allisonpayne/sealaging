#!/bin/bash

#SBATCH --job-name=MRtest
#SBATCH --mail-user=alrpayne@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2
#SBATCH --ntasks=1
#SBATCH --mem=20M
#SBATCH --cpus-per-task=1

module load miniconda3.9
conda activate sealaging
cd sealaging
Rscript cluster/markrecapture.R test 
