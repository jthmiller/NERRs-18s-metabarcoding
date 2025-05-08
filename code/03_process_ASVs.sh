### Idea: use indicator species 




############################################################################################
############################################################################################
### Filter table
### frequency and feature thresholds determined from plotting distribution of sample values
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_table.qza \
  --p-min-frequency 500 \
  --p-max-frequency 500000 \
  --p-min-features 10 \
  --m-metadata-file metadata/metadata.tsv \
  --o-filtered-table results/NERRS_18s_samples-table.qza

#### Filter table to euks only
qiime taxa filter-table \
    --i-table results/NERRS_18s_samples-table.qza \
    --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-include d__Eukaryota \
    --o-filtered-table results/NERRS_18s_euks_samples-table.qza

## Filter humans
qiime taxa filter-table \
    --i-table results/NERRS_18s_euks_samples-table.qza \
    --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-exclude "s__Homo_sapiens" \
    --o-filtered-table results/NERRS_18s_euks_hum_samples-table.qza 

# conda activate qiime2-amplicon-2024.5
## Filter high copy number taxa
qiime taxa filter-table \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
    --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-exclude p__Bacillariophyta,p__Ciliophora,c__Dinophyceae \
    --o-filtered-table results/NERRS_18s_euks_copy-Num-filtered-table.qza 


## Filter low frequency features
qiime feature-table filter-features-conditionally \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --p-prevalence 0.01 \
  --p-abundance 0.01 \
  --o-filtered-table results/NERRS_18s_euks_hum_freq-table.qza

### BARPLOT
qiime taxa barplot \
     --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
     --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
     --o-visualization results/NERRS_18s_taxa-barplot \
     --m-metadata-file metadata/metadata.tsv

### BARPLOT
qiime taxa barplot \
     --i-table results/NERRS_18s_euks_hum_freq-table.qza \
     --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
     --o-visualization results/NERRS_18s_freq-filtered_taxa-barplot \
     --m-metadata-file metadata/metadata.tsv


## OTUs
qiime taxa collapse \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 5 \
  --o-collapsed-table results/NERRS_18s_euks_hum_family-table.qza

qiime taxa collapse \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --i-taxonomy results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
  --p-level 7 \
  --o-collapsed-table results/NERRS_18s_euks_hum_genus-table.qza

qiime feature-table relative-frequency \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --o-relative-frequency-table results/relative_NERRS_18s_euks_hum_freq-table.qza

############################################################################################
############################################################################################





############################################################################################
############################################################################################
### ASV Core Metrics
#### Tree of all ASVs
qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences NERRS_18s_5_20_24_rep-seqs.qza \
   --o-alignment NERRS_18s_5_20_24_aligned-rep-seqs \
   --o-masked-alignment NERRS_18s_5_20_24_masked-aligned-rep-seqs.qza \
   --o-tree NERRS_18s_5_20_24_unrooted-tree.qza \
   --o-rooted-tree NERRS_18s_5_20_24_rooted-tree.qza \
   --p-n-threads 24

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/metadata.tsv \
    --output-dir results/core-metrics-results

qiime diversity alpha-phylogenetic \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza \
  --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity results/core-metrics-results/faith_pd

qiime diversity alpha-rarefaction \
    --i-table results/NERRS_18s_euks_hum_samples-table.qza \
    --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
    --p-max-depth 300000 \
    --m-metadata-file metadata/metadata.tsv \
    --p-min-depth 1000 \
    --p-steps 15 \
    --o-visualization results/core-metrics-results/alpha-rarefaction

qiime diversity alpha-group-significance \
    --i-alpha-diversity results/core-metrics-results/faith_pd.qza \
    --m-metadata-file metadata/metadata.tsv \
    --o-visualization results/core-metrics-results/alpha-group-significance

qiime diversity alpha-group-significance \
  --i-alpha-diversity results/core-metrics-results/observed_features_vector.qza \
  --m-metadata-file metadata/metadata.tsv \
  --o-visualization results/core-metrics-results/alpha-group-sig-obs-feats.qzv

############################################################################################
############################################################################################


# rsync -chavzP --stats gracem25@cedar.alliancecan.ca:test.fq.gz /home-wd/home/unhAW/jtmiller/watts/raw-data/cobb.sr.unh.edu/managed/230728_A01346_0115_AHM35KDRX2_16Mer072823-AW-NERRPBSF-MFNX071123/reads/PBBPtb051023_S4_L001_R2_001.fastq.gz
# rsync -avz -e ssh PBBPtb051023_S4_L001_R2_001.fastq.gz gracem25@cedar.alliancecan.ca:test.fq.gz






