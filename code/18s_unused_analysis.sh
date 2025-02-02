### 18s_unused_analysis.sh




## Here we examine how variance in Shannon diversity and other metadata changes across time (set with the state-column parameter) in the ECAM cohort, both in groups of samples (interactively selected as described below) and in individual subjects (set with the individual-id-column parameter).




## Longitudinal stats, including LME models
## https://docs.qiime2.org/2020.2/tutorials/longitudinal/
for rgn in Gulf NE N-Pacific Pacific-Island; do 



done
#### Output: rgn_subject_biplot.qzv

#Feature volatility analysis
#This pipeline identifies features that are predictive of a numeric metadata column, “state_column” (e.g., time), and plots their relative frequencies across states using interactive feature volatility plots (only important features are plotted). A supervised learning regressor is used to identify important features and assess their ability to predict sample states. state_column will typically be a measure of time, but any numeric metadata column can be used and this is not strictly a longitudinal method, unless if the individual_id_column parameter is used (in which case feature volatility plots will contain per-individual spaghetti lines, as described above).












##https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/030-tutorial-downstream/080-longitudinal.html



## state column will be on the x-axis, diversity on the y.

--p-formula 


salinity, random effect



Formulae will be in the
format "a ~ b + c", where "a" is the metric
(dependent variable) and "b" and "c" are independent
covariates. Use "+" to add a variable; "+ a:b" to
add an interaction between variables a and b; "*" to
include a variable and all interactions; and "-" to
subtract a particular term (e.g., an interaction
term).






qiime longitudinal maturity-index \
  --i-table ${rgn}_NERRS_18s-table.qza \
  --m-metadata-file qiime-corrected-sample-metadata.tsv \
  --p-state-column Quarter_num \
  --p-group-by delivery \
  --p-individual-id-column Site_Corrected \
  --p-control Vaginal \
  --p-test-size 0.4 \
  --p-stratify \
  --p-random-state 1010101 \
  --output-dir maturity

done


Non-parametric microbial interdependence test (NMIT)
Within microbial communities, microbial populations do not exist in isolation but instead form complex ecological interaction webs. Whether these interdependence networks display the same temporal characteristics within subjects from the same group may indicate divergent temporal trajectories. NMIT evaluates how interdependencies of features (e.g., microbial taxa, sequence variants, or OTUs) within a community might differ over time between sample groups. NMIT performs a nonparametric microbial interdependence test to determine longitudinal sample similarity as a function of temporal microbial composition. For each subject, NMIT computes pairwise correlations between each pair of features. Between-subject distances are then computed based on a distance norm between each subject’s microbial interdependence correlation matrix. For more details and citation, please see Zhang et al., 2017.

qiime longitudinal nmit \
  --i-table ecam-table-taxa.qza \
  --m-metadata-file ecam-sample-metadata.tsv \
  --p-individual-id-column studyid \
  --p-corr-method pearson \
  --o-distance-matrix nmit-dm.qza













qiime empress community-plot\
    --i-tree qiime2-moving-pictures-tutorial/phylo-tree.qza\
    --i-feature-table qiime2-moving-pictures-tutorial/phylo-table.qza\
    --i-pcoa qiime2-moving-pictures-tutorial/phylo-ordination.qza\
    --m-sample-metadata-file qiime2-moving-pictures-tutorial/sample-metadata.tsv\
    --m-feature-metadata-file qiime2-moving-pictures-tutorial/phylo-taxonomy.qza\
    --p-filter-missing-features\
    --p-number-of-features 50\
    --o-visualization qiime2-moving-pictures-tutorial/phylo-empress.qzv




qiime gemelli joint-rpca \
    --i-tables ${rgn}_NERRS_18s-table.qza \
    --m-sample-metadata-file qiime-corrected-sample-metadata-region.tsv \
    --p-train-test-column  'train_test' \
    --p-max-iterations 15 \
    --p-min-feature-frequency 5 \
    --o-biplot multi-omics-10333/results/joint_biplot.qza\
    --o-distance-matrix multi-omics-10333/results/joint_distance_matrix.qza\
    --o-cross-validation-error multi-omics-10333/results/cross_validation_error.qza



 qiime dev refresh-cache

for rgn in Gulf NE N-Pacific Pacific-Island; do 

     qiime deicode rpca \
    --i-table ${rgn}_NERRS_18s-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot ${rgn}_deicode-rpca-ordination.qza \
    --o-distance-matrix ${rgn}_deicode-rpca-distance.qza
  
  qiime emperor biplot \
    --i-biplot ordination.qza \
    --m-sample-metadata-file sample-metadata.tsv \
    --m-feature-metadata-file taxonomy.qza \
    --o-visualization biplot.qzv \
    --p-number-of-features 8


