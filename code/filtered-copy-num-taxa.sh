 
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=14 --time=05:00:00 --pty bash -i
conda activate /mnt/home/watts/jtm1171/.conda/envs/qiime2-amplicon-2024.5
##conda activate qiime2-amplicon-2024.5


srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=14 --time=05:00:00 --pty bash -i

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

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column salinity \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-ASVs/Salinity


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

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column salinity \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/Salinity





############### ASVs ##### 
qiime sample-classifier regress-samples \
  --i-table results/NERRS_18s_euks_SILVA_copy-Num-filtered-table.qza \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-column Sal_Min \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-ASVs/regress-Sal_Min

for file in results/copy-filtered/sample-classifier-ASVs/regress-Sal_Min/*qzv ; do cp $file ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-ASVs/Sal_Min-$(basename $file); done



qiime sample-classifier regress-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-column Sal_Min \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min

qiime sample-classifier heatmap \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --i-importance results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/feature_importance.qza \
  --m-sample-metadata-file metadata-noNA.tsv \
  --m-sample-metadata-column Sal_Min \
  --p-group-samples \
  --p-feature-count 30 \
  --o-filtered-table results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/Sal_Min-important-feature-table-top-30.qza \
  --o-heatmap results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/Sal_Min-important-feature-heatmap.qzv

qiime sample-classifier heatmap \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --i-importance results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/feature_importance.qza \
  --m-sample-metadata-file metadata-noNA.tsv \
  --m-sample-metadata-column salinity \
  --p-group-samples \
  --p-feature-count 30 \
  --o-filtered-table results/copy-filtered/sample-classifier-OTUs/regress-salinity/salinity-important-feature-table-top-30.qza \
  --o-heatmap results/copy-filtered/sample-classifier-OTUs/regress-salinity/salinity-important-feature-heatmap.qzv


for file in results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/*qzv ; do cp $file ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/Sal_Min-$(basename $file); done


here 
qiime sample-classifier regress-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-column pH_Min \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/regress-pH_Min

qiime sample-classifier heatmap \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --i-importance results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/feature_importance.qza \
  --m-sample-metadata-file metadata-noNA.tsv \
  --m-sample-metadata-column pH_Min \
  --p-group-samples \
  --p-feature-count 30 \
  --o-filtered-table results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/pH_Min-important-feature-table-top-30.qza \
  --o-heatmap results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/pH_Min-important-feature-heatmap.qzv

for file in results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/*qzv ; do cp $file ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/pH_Min-$(basename $file); done


qiime sample-classifier regress-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-column pH_Min \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/regress-pH_Min



### Do temp, nutrients 


### Other variables?

qiime sample-classifier regress-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-column pH_Min \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 20 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/regress-pH_Min





cp results/copy-filtered/sample-classifier-OTUs/regress-Sal_Min/model_summary.qzv \
~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/Sal_Min-model_summary.qzv

cp results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/accuracy_results.qzv \
~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/pH_Min-accuracy_results.qzv

for file in results/copy-filtered/sample-classifier-OTUs/regress-pH_Min/* ; do cp $file ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/pH_Min-$(basename $file); done





qiime longitudinal feature-volatility \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza \
  --m-metadata-file  metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-n-estimators 10 \
  --p-random-state 17 \
  --output-dir results/core-metrics-results/phylogenetic/longitudinal-genus
 


cp results/copy-filtered/sample-classifier-OTUs/Region/heatmap.qzv ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/Region-heatmap.qzv
cp results/copy-filtered/sample-classifier-OTUs/Site_Corrected/heatmap.qzv ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/Site-heatmap.qzv
cp results/copy-filtered/sample-classifier-OTUs/NERR/heatmap.qzv ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/NERR-heatmap.qzv



### can i run this with CTF? 

qiime gemelli ctf\
    --i-table results/NERRS_18s_euks_copy-Num-filtered-table.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/gemelli-ctf/otus


qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata-noNA.tsv \
  --m-metadata-file  results/gemelli-ctf-otu ### core-diversity/shannon_vector.qza \
  --p-metric shannon \
  --p-group-columns NERR,Region,Ocean,North_South,salinity \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --o-visualization linear-mixed-effects.qzv







qiime gemelli ctf\
    --i-table results/NERRS_18s_euks_copy-Num-filtered-table.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/gemelli-ctf/asvs


  # Step1: Generate first-distances (example distance from final)
qiime longitudinal first-distances \
  --i-distance-matrix results/gemelli-ctf/asvs/distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --p-state-column Quarter_num \
  --p-individual-id-column Site_Corrected \
  --p-replicate-handling random \
  --o-first-distances results/gemelli-ctf/asvs/first-distances.qza

# Step 2: Run LME on first-distance output
## 

python code/18s_format-subject-metadata-gemelli.py metadata.tsv
metadata-subject-metadata.tsv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file results/gemelli-ctf/asvs/first-distances.qza \
  --m-metadata-file metadata.tsv \
  --p-metric Distance \
  --p-state-column Quarter_num \
  --p-individual-id-column host_subject_id \
  --p-group-columns NERR,Region,Ocean,North_South,salinity\
  --p-formula "Distance ~  Quarter_num * salinity" \
  --o-visualization  results/gemelli-ctf/asvs/first-distances-salinity.qzv

#Step 3: Optional volatility plot of the first distances
qiime longitudinal volatility \
  --m-metadata-file metadata.tsv \
  --m-metadata-file results/gemelli-ctf/asvs/first-distances.qza \
  --p-default-metric Distance \
  --p-default-group-column ibd \
  --p-state-column timepoint \
  --p-individual-id-column host_subject_id \
  --o-visualization results/gemelli-ctf/asvs/volatility.qzv




## Filter out controls
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "NERR='SF'" \
  --o-filtered-table results/NERRS_18s_SF-table.qza


mkdir 

qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_SF-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Site_Corrected \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/Site_Corrected_SF


cp \
  results/copy-filtered/sample-classifier-OTUs/Site_Corrected_SF/*.qzv \
  ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/SF/




## Filter out controls
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_copy-Num_genus-table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "NERR='GTM'" \
  --o-filtered-table results/NERRS_18s_GTM-table.qza

  qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_GTM-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column NERR_SITE_QTR \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-missing-samples  ignore \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/SITEQTR_Corrected_GTM


mkdir ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM
cp \
  results/copy-filtered/sample-classifier-OTUs/Site_Corrected_GTM/*.qzv \
  ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM/

stratify=False



## Filter out diatoms and dinoflagellates
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "NERR='GTM'" \
  --o-filtered-table results/NERRS_18s_GTM-all-table.qza

qiime taxa collapse \
  --i-table results/NERRS_18s_GTM-all-table.qza \
  --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 7 \
  --o-collapsed-table results/NERRS_18s_GTM-g_OTU-table.qza

  qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_GTM-all-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Site_Corrected \
  --p-n-jobs 12 \
  --p-optimize-feature-selection \
  --p-missing-samples  ignore \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/SITE_Corrected_GTM_unfiltered

mkdir ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM_unfiltered
cp \
  results/copy-filtered/sample-classifier-OTUs/SITE_Corrected_GTM_unfiltered/*.qzv \
  ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM_unfiltered/


NERR_SITE_QTR


  qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_GTM-all-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column NERR_SITE_QTR \
  --p-optimize-feature-selection \
  --p-missing-samples  ignore \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/GTM_QTR_unfiltered



  qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_GTM-g_OTU-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Site_Corrected \
  --p-optimize-feature-selection \
  --p-missing-samples  ignore \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/GTM_unfiltered




  qiime sample-classifier classify-samples \
  --i-table results/NERRS_18s_GTM-g_OTU-table.qza \
  --m-metadata-file  metadata.tsv \
  --m-metadata-column Quarter_txt \
  --p-optimize-feature-selection \
  --p-missing-samples  ignore \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 123 \
  --output-dir results/copy-filtered/sample-classifier-OTUs/GTM_quarter

mkdir ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM_quarter
cp \
  results/copy-filtered/sample-classifier-OTUs/GTM_quarter/*.qzv \
  ~/jthmiller.github.io/files/results/nerrs/copy-filtered/sample-classifier-OTUs/GTM_quarter/
