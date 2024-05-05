# You should have already cloned the reo
# git clone https://github.com/FlukeAndFeather/sealaging

# Set up conda environment
module load miniconda3.9
conda create -n sealaging r-essentials r-base

# Install packages
conda activate sealaging
conda install r-posterior
conda install r-rstan

# Clone repo
echo Final step: transfer data files from your computer.
echo scp -r \<path-to-data\> $USER@hb.ucsc.edu:sealaging/data
