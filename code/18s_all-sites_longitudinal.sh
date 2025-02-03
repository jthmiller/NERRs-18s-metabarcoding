### All sites longitudinal
# qiime filter
# qiime diversity core-metrics-phylogenetic
# qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
# qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
# qiime longitudinal feature-volatility

############################################################################################################
### Filter table
qiime feature-table filter-samples \
  --i-table NERRS_18s_table.qza \
  --m-metadata-file qiime-swmp-sample-metadata.tsv \
  --o-filtered-table NERRS_18s_9_12_24_filtered-table.qza

#### Filter table to euks only
qiime taxa filter-table \
    --i-table NERRS_18s_9_12_24_filtered-table.qza \
    --i-taxonomy NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-include d__Eukaryota \
    --o-filtered-table NERRS_18s_9_12_24_euks_filtered-table.qza

## Filter humans
qiime taxa filter-table \
    --i-table NERRS_18s_9_12_24_euks_filtered-table.qza \
    --i-taxonomy NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-exclude "s__Homo_sapiens" \
    --o-filtered-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza 

qiime feature-table filter-features-conditionally \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --p-prevalence 0.01 \
  --p-abundance 0.01 \
  --o-filtered-table NERRS_18s_9_12_24_euks_hum_freq-table.qza

qiime taxa collapse \
  --i-table NERRS_18s_9_12_24_euks_hum_freq-table.qza \
  --i-taxonomy NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 5 \
  --o-collapsed-table NERRS_18s_euks_hum_family-table.qza

qiime taxa collapse \
  --i-table NERRS_18s_9_12_24_euks_hum_freq-table.qza \
  --i-taxonomy NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 7 \
  --o-collapsed-table NERRS_18s_euks_hum_genus-table.qza

qiime feature-table relative-frequency \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --o-relative-frequency-table relative_NERRS_18s_9_12_24_euks_hum_freq-table.qza


#qiime feature-table relative-frequency \
#  --i-table NERRS_18s_9_12_24_euks_hum_freq-family-table.qza \
#  --o-relative-frequency-table relative_NERRS_18s_9_12_24_euks_hum_freq-family-table.qza
############################################################################################################



############################################################################################################
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny qiime-files/input/NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table qiime-files/input/NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata.tsv \
    --output-dir qiime-files/core-diversity-phylogenetic

############################################################################################################

############################################################################################################

  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table qiime-files/input/NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
      --i-phylogeny qiime-files/input/NERRS_18s_9_12_24_rooted-tree.qza \
      --m-taxonomy-file qiime-files/input/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 50 \
      --p-min-sample-count 500 \
      --o-biplot qiime-files/core-diversity-phylogenetic/phylo-ordination.qza \
      --o-distance-matrix qiime-files/core-diversity-phylogenetic/phylo-distance.qza \
      --o-counts-by-node-tree qiime-files/core-diversity-phylogenetic/phylo-tree.qza \
      --o-counts-by-node qiime-files/core-diversity-phylogenetic/phylo-table.qza \
      --o-t2t-taxonomy qiime-files/core-diversity-phylogenetic/phylo-taxonomy.qza

  qiime empress community-plot\
    --i-tree qiime-files/core-diversity-phylogenetic/phylo-tree.qza\
    --i-feature-table qiime-files/core-diversity-phylogenetic/phylo-table.qza\
    --i-pcoa qiime-files/core-diversity-phylogenetic/phylo-ordination.qza\
    --m-sample-metadata-file qiime-swmp-corrected-sample-metadata.tsv\
    --m-feature-metadata-file qiime-files/core-diversity-phylogenetic/phylo-taxonomy.qza\
    --p-filter-missing-features\
    --p-number-of-features 50\
    --o-visualization qiime-files/core-diversity-phylogenetic/phylo-empress.qzv


### Rerun this
qiime diversity beta-group-significance \
    --i-distance-matrix qiime-files/core-diversity-phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column salinity \
    --p-method permanova \
    --o-visualization qiime-files/core-diversity-phylogenetic/phylo-salinity_significance.qzv
   
## Qu
## qiime 2020.8 needed

# qiime qurro loading-plot \
#     --i-ranks core-diversity-phylogenetic/phylo-ordination.qza \
#     --i-table core-diversity-phylogenetic/phylo-table.qza \
#     --m-sample-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
#     --m-feature-metadata-file core-diversity-phylogenetic/phylo-taxonomy.qza \
#     --o-visualization core-diversity-phylogenetic/phylo-qurro_plot.qzv   
############################################################################################################



conda create -n qurro python qurro qiime2

############################################################################################################
### Feature volitility on entire dataset without Bacteria and Humans
### currently doesnt work on filtered data

qiime longitudinal feature-volatility \
  --i-table qiime-files/input/NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir qiime-files/core-diversity-phylogenetic/longitudinal-filtered

qiime longitudinal feature-volatility \
  --i-table qiime-files/input/NERRS_18s_9_12_24_filtered-table.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir qiime-files/core-diversity-phylogenetic/longitudinal

### Feature volitility on entire dataset without Bacteria and Humans low freq features removed
### This doesnt factor in the relationship between sites within a NERR
qiime longitudinal feature-volatility \
  --i-table qiime-files/input/NERRS_18s_euks_hum_genus-table.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir qiime-files/core-diversity-phylogenetic/longitudinal-genus


  
############################################################################################################


############################################################################################################
qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
### gemelli ctf
qiime gemelli ctf\
    --i-table qiime-files/input/NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file qiime-files/input/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir core-diversity-phylogenetic/gemelli-ctf-asv

qiime longitudinal volatility \
    --i-table qiime-files/input/relative_NERRS_18s_9_12_24_euks_hum_freq-table.qza \
    --p-state-column Quarter_num \
    --m-metadata-file core-diversity-phylogenetic/gemelli-ctf-asv/state_subject_ordination.qza \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization core-diversity-phylogenetic/gemelli-ctf-asv/rf-state_subject_ordination.qzv      























qiime gemelli ctf\
    --i-table NERRS_18s_euks_hum_genus-table.qza \
    --m-sample-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
    --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir core-diversity-phylogenetic/gemelli-ctf-genus





qiime longitudinal volatility \
    --i-table family-rf-table.qza \
    --p-state-column Quarter_num \
    --m-metadata-file core-diversity-phylogenetic/gemelli-ctf/state_subject_ordination.qza \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization core-diversity-phylogenetic/gemelli-ctf/rf-state_subject_ordination.qzv      

qiime longitudinal volatility \
  --i-table NERRS_18s_euks_hum_genus-table.qza \
  --p-state-column Quarter_num \
  --m-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
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
  --m-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
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
  --m-metadata-file metadata.tsv \
  --m-metadata-file qiime-files/all-sites/core-diversity-phylogenetic/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-columns NERR,Region,Ocean,North_South,salinity \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --o-visualization qiime-files/visualization-files/all-sites/linear-mixed-effects.qzv



qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.tsv \
  --m-metadata-file qiime-files/all-sites/core-diversity-phylogenetic/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-columns region \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --o-visualization qiime-files/visualization-files/all-sites/linear-mixed-effects-region.qzv



qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.tsv \
  --m-metadata-file qiime-files/all-sites/core-diversity-phylogenetic/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-columns salinity \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --o-visualization qiime-files/visualization-files/all-sites/linear-mixed-effects-salinity.qzv


"/home/users/jtm1171/.conda/envs/qiime2-amplicon-2024.5/lib/python3.9/site-packages/q2_types/sample_data/_transformer.py:27: FutureWarning: errors='ignore' is deprecated and will raise in a future version. Use to_numeric without passing `errors` and catch exceptions explicitly instead
  df[cols] = df[cols].apply(pd.to_numeric, errors='ignore')"


### Rerun this
qiime diversity beta-group-significance \
    --i-distance-matrix qiime-files/all-sites/core-diversity-phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column salinity \
    --p-method permanova \
    --o-visualization qiime-files/visualization-files/all-sites/phylo-salinity_significance.qzv


qiime diversity core-metrics \
    --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
    --output-dir core-diversity/

