# NERRs-18s-metabarcoding
Metabarcoding analysis for NERRs sites sampled quarterly. Metadata for all sites [here](metadata.tsv), including SWMP data when available. Downloaded from https://cdmo.baruch.sc.edu
## OTU Abundance 

![gulf](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/gulf-barplots-sample.png?raw=true)


## Core diversity metrics 
PCoA is performed on distance matrices for the metrics below (seems to better handle missing data than PCA does). 

More complete descriptions here: https://docs.onecodex.com/en/articles/4150649-beta-diversity
* Alpha diversity
    * Shannon’s diversity index (a quantitative measure of community richness)
    * Observed OTUs (a qualitative measure of community richness)
    * Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic relationships between the features)
    * Evenness (or Pielou’s Evenness; a measure of community evenness)
* Beta diversity
    * Jaccard distance (a qualitative measure of community dissimilarity. Qualitative - presence / absence - percentage of taxa not found in both samples) [All Samples](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-diversity-phylogenetic/jaccard_emperor.qzv)
    * Bray-Curtis distance (a quantitative measure of community dissimilarity. Takes into consideration abundance and presence absence) [All Samples](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-diversity-phylogenetic/bray_curtis_emperor.qzv)

    * Unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features. Percentage of phylogenetic branch length not found in both samples)
    [All Samples](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-diversity-phylogenetic/unweighted_unifrac_emperor.qzv)
    * Weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features. Similar to Bray-Curtis but takes into consideration phylogenetic relationships)
    [All Samples](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-diversity-phylogenetic/weighted_unifrac_emperor.qzv)



## These files still need to be updated with new salinity annotation 
[phylo-empress](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/phylo-empress.qzv)

[linear-mixed-effects-by-region](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/linear-mixed-effects-region.qzv)

[linear-mixed-effects-by-salinity](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/linear-mixed-effects-salinity.qzv)

[phylo-salinity_significance](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/phylo-salinity_significance.qzv)

[rf-state_subject_ordination](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/rf-state_subject_ordination.qzv)

[volatility_plot](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/volatility_plot.qzv)

[genus-level_volatility](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/genus-level_volatility_plot.qzv)

[ASV-level_volatility](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/all-sites/asv-level_volatility_plot.qzv)

### New England NERRS results
[NE_subject_biplot](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/regional/NE_with-repl/NE_gemelli-ctf/NE_subject_biplot.qzv)

[state_subject_ordination](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/regional/NE_with-repl/NE_gemelli-ctf/state_subject_ordination.qzv)

[accuracy_results](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/regional/NE_with-repl/NE_ecam-feat-volatility/accuracy_results.qzv)



### Unifrac PCoA performed on Unweighted UniFrac distance matrix 
![unifrac](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/unifrac_salinity_all-sites.png?raw=true)
Samples colored by minimum salinity from SWMP collected data within X days of eDNA sample collection. This is an interactive plot that can be found [here](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/core-diversity-phylogenetic/weighted_unifrac_emperor.qzv)


### Results:  

Bar-plot images of sites broken up by regions [here](images/barplots)  

Network plots of sites broken up by regions [here](images/network-plots/)  

We also perfomed PERMANOVA to test for differences between groups.

## Longitudinal Analysis
From the tutorial: Repeat measure experimental designs (e.g. time series) are a valid and powerful method to control for inter-individual variation. However, conventional dimensionality reduction methods can not account for the high-correlation of each subject to itself at a later time point. 

Longitudinal analysis can occur on features (ASVs), taxonomic assignment, or on phylogenetic clusters. More info and tutorials for longitudinal stats, including LME models can be found [here](https://docs.qiime2.org/2020.2/tutorials/longitudinal/)

#### Compositional tensor factorization (CTF) 
In order to account for the correlation among samples from the same subject we will employ compositional tensor factorization (CTF). CTF builds on the ability to account for compositionality and sparsity using the robust center log-ratio transform covered in the RPCA tutorial (found here) but restructures and factors the data as a tensor. Here we will run CTF through gemelli and explore/interpret the different results.

The package 'gemelli' for qiime2 is used to perform CTF.  The gemelli tutorial can be found [here](https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/IBD-Tutorial-QIIME2-CLI.md) and at doi: 10.1038/s41587-020-0660-7


[link to phylo-salinity_significance](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/all-sites/phylo-salinity_significance.qzv)



- Treat each NERR_SITE as a subject
- Treat each NERR_SITE_QTR as a time

The volatility visualizer generates interactive line plots that allow us to assess how volatile a dependent variable is over a continuous, independent variable (e.g., time) in one or more groups. 

Multiple metadata files (including alpha and beta diversity artifacts) and FeatureTable[RelativeFrequency] tables can be used as input, and in the interactive visualization we can select different dependent variables to plot on the y-axis.




#### state_subject_ordination

The y-axis in the subject trajectory is a PC axis like a conventional ordination (i.e. PCoA) and the x-axis is time.
The interpretation is also similar to a conventional ordination scatter plot -- where the larger the distance is between subjects at each time point the greater the difference in their communities. 

Here we can see that CTF can effectively show a difference between sites subjects across time.

![ctf-volitility](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/NE-state-subject-ordination-ctf-sample.png?raw=true)

#### subject_biplot.qzv

We can also see that the IBD grouping is separated entirely along the first PC (axis 1). We can now use Qurro to explore the feature loading partitions (arrows) in this biplot as a log-ratio of the original table counts. This allows us to relate these low-dimensional representations back to our original data. Additionally, log-ratios provide a nice set of data points for additional analysis such as LME models.

#### Feature volatility  


![fv](https://github.com/jthmiller/NERRs-18s-metabarcoding/blob/main/images/sample-plots/volatility-control-chart.png?raw=true)



