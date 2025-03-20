

## ## NERRs 18s Final Analysis
conda config --prepend pkgs_dirs /home/users/jtm1171/.conda/pkgs
conda config --prepend envs_dirs /home/users/jtm1171/.conda/envs


####### NEW ENVIRONMENT FOR diversity metrics gemelli
conda env create -n qiime2-div --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
conda activate qiime2-div
conda install bioconda::bioconductor-phyloseq r-qiime2r conda-forge::r-tidyverse conda-forge::r-devtools
pip install empress
qiime dev refresh-cache
## qurro breaks conda env. Dont install


conda list --revisions
conda install --revision 0



## ##### NEW ENVIRONMENT FOR Qiime2R
 conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
 conda activate qiime2-amplicon-2024.5
 conda install bioconda::bioconductor-phyloseq r-qiime2r conda-forge::r-tidyverse conda-forge::r-devtools
 
 
 conda install r-qiime2r
 conda install bioconda::bioconductor-phyloseq
 conda install conda-forge::r-tidyverse
 conda install conda-forge::r-devtools