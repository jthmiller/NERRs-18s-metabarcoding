### Subset 18 NERRs data to perform PCoAs by region

cd /home/users/jtm1171/NERRs/18s
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

## Filter by region
qiime feature-table filter-samples \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --m-metadata-file qiime-corrected-sample-metadata-region.tsv \
  --p-where "NERR='HE'" \
  --o-filtered-table Pacific-Island_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --m-metadata-file qiime-corrected-sample-metadata-region.tsv \
  --p-where "Region='Gulf'" \
  --o-filtered-table Gulf_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --m-metadata-file qiime-corrected-sample-metadata-region.tsv \
  --p-where "Region='NE'" \
  --o-filtered-table NE_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table NERRS_18s_9_12_24_euks_hum_filtered-table.qza \
  --m-metadata-file qiime-corrected-sample-metadata-region.tsv \
  --p-where "Region='N-Pacific'" \
  --o-filtered-table N-Pacific_NERRS_18s-table.qza
