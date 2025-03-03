### All sites longitudinal
# qiime filter
# qiime diversity core-metrics-phylogenetic
# qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
# qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
# qiime longitudinal feature-volatility

############################################################################################################
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata/metadata.tsv \
    --output-dir results/core-metrics-results/phylogenetic/

############################################################################################################

############################################################################################################

  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table results/NERRS_18s_euks_hum_samples-table.qza \
      --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
      --m-taxonomy-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 50 \
      --p-min-sample-count 500 \
      --o-biplot results/core-metrics-results/phylogenetic/phylo-ordination.qza \
      --o-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
      --o-counts-by-node-tree results/core-metrics-results/phylogenetic/phylo-tree.qza \
      --o-counts-by-node results/core-metrics-results/phylogenetic/phylo-table.qza \
      --o-t2t-taxonomy results/core-metrics-results/phylogenetic/phylo-taxonomy.qza

  qiime empress community-plot\
    --i-tree results/core-metrics-results/phylogenetic/phylo-tree.qza\
    --i-feature-table results/core-metrics-results/phylogenetic/phylo-table.qza\
    --i-pcoa results/core-metrics-results/phylogenetic/phylo-ordination.qza\
    --m-sample-metadata-file metadata/metadata.tsv\
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza\
    --p-filter-missing-features\
    --p-number-of-features 50\
    --o-visualization results/core-metrics-results/phylogenetic/phylo-empress.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata/metadata.tsv \
    --m-metadata-column salinity \
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-salinity_significance.qzv
   
## Qu
## qiime 2020.8 needed

qiime qurro loading-plot \
    --i-ranks results/core-metrics-results/phylogenetic/phylo-ordination.qza \
    --i-table results/core-metrics-results/phylogenetic/phylo-table.qza \
    --m-sample-metadata-file metadata/metadata.tsv \
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-qurro_plot.qzv   
############################################################################################################



conda create -n qurro python qurro qiime2

############################################################################################################
### Feature volitility on entire dataset without Bacteria and Humans
### currently doesnt work on filtered data

qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --m-metadata-file metadata/metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-filtered

qiime longitudinal feature-volatility \
  --i-table NERRS_18s_9_12_24_filtered-table.qza \
  --m-metadata-file metadata/metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal

### Feature volitility on entire dataset without Bacteria and Humans low freq features removed
### This doesnt factor in the relationship between sites within a NERR
qiime longitudinal feature-volatility \
  --i-table NERRS_18s_euks_hum_genus-table.qza \
  --m-metadata-file metadata/metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-genus


  
############################################################################################################


############################################################################################################
qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
### gemelli ctf
qiime gemelli ctf\
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --m-sample-metadata-file metadata/metadata.tsv \
    --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/core-metrics-results/phylogenetic/gemelli-ctf-asv

qiime longitudinal volatility \
    --i-table relative_NERRS_18s_9_12_24_euks_hum_freq-table.qza \
    --p-state-column Quarter_num \
    --m-metadata-file results/core-metrics-results/phylogenetic/gemelli-ctf-asv/state_subject_ordination.qza \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf-asv/rf-state_subject_ordination.qzv      























qiime gemelli ctf\
    --i-table NERRS_18s_euks_hum_genus-table.qza \
    --m-sample-metadata-file metadata/metadata.tsv \
    --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/core-metrics-results/phylogenetic/gemelli-ctf-genus





qiime longitudinal volatility \
    --i-table family-rf-table.qza \
    --p-state-column Quarter_num \
    --m-metadata-file results/core-metrics-results/phylogenetic/gemelli-ctf/state_subject_ordination.qza \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf/rf-state_subject_ordination.qzv      

qiime longitudinal volatility \
  --i-table NERRS_18s_euks_hum_genus-table.qza \
  --p-state-column Quarter_num \
  --m-metadata-file metadata/metadata.tsv \
  --p-individual-id-column subject_id \
  --p-default-group-column NERR \
  --o-visualization genus-all-volatility-plot-1.qzv



#### FEATURE VOLITILITY
qiime longitudinal feature-volatility \
  --i-table filtered-family-table.qza \
  --m-metadata-file qiime-swmp-corrected-subject-metadata.tsv  \
  --p-state-column Quarter_num \
  --p-individual-id-column subject_id \
  --output-dir family-longitudinal-feature-volatility

qiime longitudinal feature-volatility \
  --i-table filtered-family-table.qza \
  --m-metadata-file metadata/metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir all_ecam-feat-volatility-family

qiime longitudinal feature-volatility \
  --i-table genus-rf-table.qza \
  --m-metadata-file sample-metadata.tsv uu-umap.qza diversity-core-metrics-phylogenetic/faith_pd_vector.qza diversity-core-metrics-phylogenetic/evenness_vector.qza diversity-core-metrics-phylogenetic/shannon_vector.qza \
  --p-state-column week-relative-to-fmt \
  --p-individual-id-column subject_id \
  --output-dir longitudinal-feature-volatility-2






## LME models 
Linear mixed effects (LME) models test the relationship between a single response variable and one or more independent variables, where observations are made across dependent samples, e.g., in repeated-measures sampling experiments. This implementation takes at least one numeric state-column (e.g., Time) and one or more comma-separated group-columns (which may be categorical or numeric metadata columns; these are the fixed effects) as independent variables in a LME model, and plots regression plots of the response variable (“metric”) as a function of the state column and each group column. Additionally, the individual-id-column parameter should be a metadata column that indicates the individual subject/site that was sampled repeatedly. 
- NERR_SITE is repeatedly measured over quarters





qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata/metadata.tsv \
  --m-metadata-file core-diversity/shannon_vector.qza \
  --p-metric shannon \
  --p-group-columns NERR,Region,Ocean,North_South,salinity \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --o-visualization linear-mixed-effects.qzv


qiime diversity core-metrics \
    --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata/metadata.tsv \
    --output-dir core-diversity/

