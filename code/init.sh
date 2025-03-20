
# You may need to install conda into your home directory before this will work. 
conda config --prepend pkgs_dirs /home/users/<your-username>/.conda/pkgs
conda config --prepend envs_dirs /home/users/<your-username>/.conda/envs

# New FOR Qiime2
conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
# you will need to confirm

conda activate qiime2-amplicon-2024.5
conda install bioconda::bioconductor-phyloseq r-qiime2r conda-forge::r-tidyverse conda-forge::r-devtools
 