 
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=14 --time=05:00:00 --pty bash -i
conda activate /mnt/home/watts/jtm1171/.conda/envs/qiime2-amplicon-2024.5
##conda activate qiime2-amplicon-2024.5

# conda activate qiime2-amplicon-2024.5
## Filter high copy number taxa
qiime taxa filter-table \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
    --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-exclude p__Bacillariophyta,p__Ciliophora,c__Dinophyceae \
    --o-filtered-table results/NERRS_18s_euks_copy-Num-filtered-table.qza 

qiime taxa filter-table \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
    --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-exclude p__Bacillariophyta,p__Ciliophora,p__Diatomea,c__Dinophyceae \
    --o-filtered-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza 

qiime taxa collapse \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 7 \
  --o-collapsed-table results/NERRS_18s_euks_copy-Num_genus-table.qza



############################################################################################################
rm -fR results/copy-filtered/core-metrics-results

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata.tsv \
    --output-dir results/copy-filtered/core-metrics-OTUs

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 1000 \
    --m-metadata-file metadata.tsv \
    --output-dir results/copy-filtered/core-metrics-ASVs


cp -R results/copy-filtered /mnt/gpfs01/home/watts/jtm1171/jthmiller.github.io/files/results/nerrs/
############################################################################################################


############################################################################################################
  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
      --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
      --m-taxonomy-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 50 \
      --p-min-sample-count 500 \
      --o-biplot results/copy-filtered/core-metrics-results/phylogenetic/phylo-ordination.qza \
      --o-distance-matrix results/copy-filtered/core-metrics-results/phylogenetic/phylo-distance.qza \
      --o-counts-by-node-tree results/copy-filtered/core-metrics-results/phylogenetic/phylo-tree.qza \
      --o-counts-by-node results/copy-filtered/core-metrics-results/phylogenetic/phylo-table.qza \
      --o-t2t-taxonomy results/copy-filtered/core-metrics-results/phylogenetic/phylo-taxonomy.qza

qiime qurro loading-plot \
    --i-ranks results/copy-filtered/core-metrics-results/phylogenetic/phylo-ordination.qza \
    --i-table results/copy-filtered/core-metrics-results/phylogenetic/phylo-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/copy-filtered/core-metrics-results/phylogenetic/phylo-taxonomy.qza \
    --o-visualization results/copy-filtered/core-metrics-results/phylogenetic/phylo-qurro_plot.qzv

  qiime empress community-plot\
    --i-tree results/copy-filtered/core-metrics-results/phylogenetic/phylo-tree.qza\
    --i-feature-table results/copy-filtered/core-metrics-results/phylogenetic/phylo-table.qza\
    --i-pcoa results/copy-filtered/core-metrics-results/phylogenetic/phylo-ordination.qza\
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/copy-filtered/core-metrics-results/phylogenetic/phylo-taxonomy.qza\
    --p-filter-missing-features \
    --p-number-of-features 50 \
    --o-visualization results/copy-filtered/core-metrics-results/phylogenetic/phylo-empress.qzv


qiime diversity beta-group-significance \
    --i-distance-matrix results/copy-filtered/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column salinity \
    --p-method permanova \
    --o-visualization results/copy-filtered/core-metrics-results/phylogenetic/phylo-salinity_significance.qzv
   
qiime diversity beta-group-significance \
    --i-distance-matrix results/copy-filtered/core-metrics-results/phylogenetic/phylo-distance.qza \
    --m-metadata-file metadata.tsv \
    --m-metadata-column  Site_Corrected\
    --p-method permanova \
    --o-visualization results/copy-filtered/core-metrics-results/phylogenetic/phylo-Site_significance.qzv

qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/copy-filtered/core-metrics-results/phylogenetic/longitudinal-filtered




###### ######


qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Region \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-ASVs/Region

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Site_Corrected \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-ASVs/Site_Corrected

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column NERR \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-ASVs/NERR


### OTUs
rm results/NERRS_18s_euks_copy-Num_genus-table.qza
qiime taxa collapse \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 7 \
  --o-collapsed-table results/NERRS_18s_euks_copy-Num_genus-table.qza

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Region \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/Region

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Site_Corrected \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/Site_Corrected

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column NERR \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/NERR










qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-genus
 