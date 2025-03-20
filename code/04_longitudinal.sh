### All sites longitudinal
# qiime filter
# qiime diversity core-metrics-phylogenetic
# qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
# qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
# qiime longitudinal feature-volatility
?? qiime longitudinal pairwise-differences



conda activate qiime2-amplicon-2024.5
############################################################################################################
# rm -fR results/core-metrics-results/core
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/core-metrics-results/NERRS_18s_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata.tsv \
    --output-dir results/core-metrics-results/core
cp results/core-metrics-results/core/*qzv /mnt/gpfs01/home/watts/jtm1171/jthmiller.github.io/files/results/nerrs/core-metrics-results/
############################################################################################################


############################################################################################################

  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table results/NERRS_18s_euks_hum_samples-table.qza \
      --i-phylogeny results/core-metrics-results/NERRS_18s_rooted-tree.qza \
      --m-taxonomy-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 50 \
      --p-min-sample-count 500 \
      --o-biplot results/core-metrics-results/phylogenetic/phylo-ordination.qza \
      --o-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
      --o-counts-by-node-tree results/core-metrics-results/phylogenetic/phylo-tree.qza \
      --o-counts-by-node results/core-metrics-results/phylogenetic/phylo-table.qza \
      --o-t2t-taxonomy results/core-metrics-results/phylogenetic/phylo-taxonomy.qza

conda activate qiime2-amplicon-2024.5
qiime qurro loading-plot \
    --i-ranks results/core-metrics-results/phylogenetic/phylo-ordination.qza \
    --i-table results/core-metrics-results/phylogenetic/phylo-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-qurro_plot.qzv

  qiime empress community-plot\
    --i-tree results/core-metrics-results/phylogenetic/phylo-tree.qza\
    --i-feature-table results/core-metrics-results/phylogenetic/phylo-table.qza\
    --i-pcoa results/core-metrics-results/phylogenetic/phylo-ordination.qza\
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza\
    --p-filter-missing-features \
    --p-number-of-features 50 \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-empress.qzv


for 
  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table results/NERRS_18s_euks_hum_samples-table.qza \
      --i-phylogeny results/core-metrics-results/NERRS_18s_rooted-tree.qza \
      --m-taxonomy-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 50 \
      --p-min-sample-count 500 \
      --o-biplot results/core-metrics-results/phylogenetic/phylo-ordination.qza \
      --o-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
      --o-counts-by-node-tree results/core-metrics-results/phylogenetic/phylo-tree.qza \
      --o-counts-by-node results/core-metrics-results/phylogenetic/phylo-table.qza \
      --o-t2t-taxonomy results/core-metrics-results/phylogenetic/phylo-taxonomy.qza

conda activate qiime2-amplicon-2024.5
qiime qurro loading-plot \
    --i-ranks results/core-metrics-results/phylogenetic/phylo-ordination.qza \
    --i-table results/core-metrics-results/phylogenetic/phylo-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-qurro_plot.qzv

  qiime empress community-plot\
    --i-tree results/core-metrics-results/phylogenetic/phylo-tree.qza\
    --i-feature-table results/core-metrics-results/phylogenetic/phylo-table.qza\
    --i-pcoa results/core-metrics-results/phylogenetic/phylo-ordination.qza\
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza\
    --p-filter-missing-features \
    --p-number-of-features 50 \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-empress.qzv


############################################################################################################
qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column salinity \
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-salinity_significance.qzv
   
qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column NERR \
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-NERR_significance.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column Region \
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-region_significance.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column Quarter_txt \
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-Quarter_significance.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column  Site_Corrected\
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-Site_significance.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix results/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column  Site_Corrected\
    --p-method permanova \
    --o-visualization results/core-metrics-results/phylogenetic/phylo-Site_significance.qzv

## mkdir  ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic 
cp results/core-metrics-results/phylogenetic/phylo* ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/
## Qu
## qiime 2020.8 needed

## qiime qurro loading-plot \
##     --i-ranks results/core-metrics-results/phylogenetic/phylo-ordination.qza \
##     --i-table results/core-metrics-results/phylogenetic/phylo-table.qza \
##     --m-sample-metadata-file metadata/metadata.tsv \
##     --m-feature-metadata-file results/core-metrics-results/phylogenetic/phylo-taxonomy.qza \
##     --o-visualization results/core-metrics-results/phylogenetic/phylo-qurro_plot.qzv   
############################################################################################################
## conda create -n qurro python qurro qiime2

############################################################################################################
### Feature volitility on entire dataset without Bacteria and Humans
rm -fR results/core-metrics-results/phylogenetic/longitudinal-filtered
rm -fR results/core-metrics-results/phylogenetic/longitudinal-genus
rm -fR results/core-metrics-results/phylogenetic/longitudinal-family

qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-filtered

### Feature volitility on entire dataset without Bacteria and Humans low freq features removed
### This doesnt factor in the relationship between sites within a NERR
qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-genus
 
qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_family-table.qza \
  --m-metadata-file  metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-family



qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file  metadata-noNA.tsv \
  --p-state-column Sal_Min \
  --p-individual-id-column NERR_SITE_QTR \
  --p-n-estimators 100 \
  --p-random-state 1 \
  --p-n-jobs 12 \
  --o-filtered-table results/ctf/NERRs_ctf_Sal-Min_filtered-table.qzv \
  --o-feature-importance results/ctf/NERRs_ctf_Sal-Min_feature-importance.qzv \
  --o-volatility-plot results/ctf/NERRs_ctf_Sal-Min_volatility-plot.qzv \
  --o-accuracy-results results/ctf/NERRs_ctf_Sal-Min_accuracy-results.qzv \
  --o-sample-estimator results/ctf/NERRs_ctf_Sal-Min_sample-estimator.qzv &

Usage: qiime diversity adonis [OPTIONS]


 cp -r results/core-metrics-results/phylogenetic/longitudinal-filtered ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/
 cp -r results/core-metrics-results/phylogenetic/longitudinal-genus ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/
 cp -r results/core-metrics-results/phylogenetic/longitudinal-family ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/


### Regional feature volitility




############################################################################################################
### to compare two specific states (quarters), use 
mkdir results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis


qiime longitudinal pairwise-distances \
  --i-distance-matrix results/core-metrics-results/core/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column Region \
  --p-state-column Quarter_num \
  --p-state-1 1 \
  --p-state-2 3 \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-visualization results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/Region_1-3_pairwise-distances.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix results/core-metrics-results/core/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column Region \
  --p-state-column Quarter_num \
  --p-state-1 2 \
  --p-state-2 4 \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-visualization results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/Region_2-4_pairwise-distances.qzv



qiime longitudinal pairwise-distances \
  --i-distance-matrix results/core-metrics-results/core/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column North_South \
  --p-state-column Quarter_num \
  --p-state-1 1 \
  --p-state-2 3 \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-visualization results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/North_South-1_3-pairwise-distances.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix results/core-metrics-results/core/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column North_South \
  --p-state-column Quarter_num \
  --p-state-1 2 \
  --p-state-2 4 \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-visualization results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/North_South-2_4-pairwise-distances.qzv


cp results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/* ~/jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/
############################################################################################################

## qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
### gemelli ctf
python code/18s_format-subject-metadata-gemelli.py metadata.tsv

qiime gemelli ctf\
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/core-metrics-results/phylogenetic/gemelli-ctf-asv

qiime qurro loading-plot\
    --i-table results/NERRS_18s_euks_hum_samples-table.qza\
    --i-ranks results/core-metrics-results/phylogenetic/gemelli-ctf-asv/subject_biplot.qza\
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza\
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf-asv/qurro.qzv




qiime emperor biplot\
    --i-biplot results/core-metrics-results/phylogenetic/gemelli-ctf-asv/subject_biplot.qza \
    --m-sample-metadata-file metadata-subject-metadata.tsv \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-number-of-features 50 \
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf-asv/gemelli-ctf-asv_subject_biplot.qzv



qiime longitudinal volatility \
    --p-state-column Quarter_num \
    --m-metadata-file results/core-metrics-results/phylogenetic/gemelli-ctf-asv/state_subject_ordination.qza \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf-asv/rf-state_subject_ordination.qzv      

qiime longitudinal volatility \
    --p-state-column Quarter_num \
    --m-metadata-file results/core-metrics-results/core/unweighted_unifrac_emperor.qzv \
    --p-individual-id-column subject_id \
    --p-default-group-column NERR \
    --p-default-metric PC1 \
    --o-visualization results/core-metrics-results/phylogenetic/gemelli-ctf-asv/rf-state_subject_ordination.qzv      


mkdir ~/jthmiller.github.io/files/results/nerrs/longitudinal 
cp results/core-metrics-results/phylogenetic/gemelli-ctf-asv/rf-state_subject_ordination.qzv ~/jthmiller.github.io/files/results/nerrs/longitudinal/


qiime longitudinal pairwise-distances \
  --i-distance-matrix results/core-metrics-results/phylogenetic/gemelli-ctf-asv/distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column North_South \
  --p-state-column Quarter_num \
  --p-state-1 2 \
  --p-state-2 4 \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-visualization results/core-metrics-results/phylogenetic/longitudinal-pairwise-dis/North_South-2_4-gemelli_ctf_pairwise-distances.qzv

results/core-metrics-results/phylogenetic/gemelli-ctf-asv/




















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

