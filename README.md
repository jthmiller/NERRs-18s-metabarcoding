# NERRs-18s-metabarcoding
Metabarcoding analysis for NERRs sites sampled quarterly 

Metadata for all sites [here](metadata.tsv)

## OTU Abundance 

![gulf](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/gulf-barplots-sample.png?raw=true)


## Core diversity metrics 
More complete descriptions here: https://docs.onecodex.com/en/articles/4150649-beta-diversity
* Alpha diversity
    * Shannon’s diversity index (a quantitative measure of community richness)
    * Observed OTUs (a qualitative measure of community richness)
    * Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic relationships between the features)
    * Evenness (or Pielou’s Evenness; a measure of community evenness)
* Beta diversity
    * Jaccard distance (a qualitative measure of community dissimilarity. Qualitative - presence / absence - percentage of taxa not found in both samples)
    * Bray-Curtis distance (a quantitative measure of community dissimilarity. Takes into consideration abundance and presence absence)
    * Unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features. Percentage of phylogenetic branch length not found in both samples)
    * Weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features. Similar to Bray-Curtis but takes into consideration phylogenetic relationships)
Ordination
- PCA is performed on the underlying data matrix
- PCoA is performed on distance matrices above (seems to better handle missing data than PCA does)
- PERMANOVA to test for differences between groups



* DEICODE is a tool box for running Robust Aitchison PCA on sparse compositional omics datasets, linking specific features to beta-diversity ordination.


![unifrac](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/unifrac_salinity_all-sites.png?raw=true)
Unifrac of all samples colored by salinity designation from SWMP collected data


The folders containing:  

Bar-plot images of sites broken up by regions [here](images/barplots)  

Network plots of sites broken up by regions [here](images/network-plots/)  


### Longitudinal Analysis
From the tutorial: Repeat measure experimental designs (e.g. time series) are a valid and powerful method to control for inter-individual variation. However, conventional dimensionality reduction methods can not account for the high-correlation of each subject to itself at a later time point. 

Longitudinal analysis can occur on features (ASVs), taxonomic assignment, or on phylogenetic clusters. 

Longitudinal stats, including LME models can be found [here](https://docs.qiime2.org/2020.2/tutorials/longitudinal/)

### Compositional tensor factorization (CTF)  doi: 10.1038/s41587-020-0660-7
The package 'gemelli' for qiime2 is used to perform CTF.  The gemelli tutorial can be found [here](https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/IBD-Tutorial-QIIME2-CLI.md)

- Treat each NERR_SITE as a subject
- Treat each NERR_SITE_QTR as a time

The volatility visualizer generates interactive line plots that allow us to assess how volatile a dependent variable is over a continuous, independent variable (e.g., time) in one or more groups. 

Multiple metadata files (including alpha and beta diversity artifacts) and FeatureTable[RelativeFrequency] tables can be used as input, and in the interactive visualization we can select different dependent variables to plot on the y-axis.





