####### Classify rep seqs ##############################################################
########################################################################################

## Try CRUX database
## taxonomy
maxaccepts=10
query_cov=0.8 
perc_identity=0.90 
weak_id=0.80
refreads=${refreads:-ref-dbs/CRUX/18S_.fasta.qza}
reftax=${reftax:-ref-dbs/CRUX/18S_taxonomy.qza}
sklearn=${sklearn:-ref-dbs/CRUX/18s_crux_nb-classifier.qza}

qiime feature-classifier classify-consensus-vsearch \
  --i-query NERRS_18s_5_20_24_rep-seqs.qza \
    --i-reference-reads ${refreads} \
    --i-reference-taxonomy  ${reftax} \
    --p-threads 41 \
    --p-perc-identity ${perc_identity} \
    --p-query-cov ${query_cov} \
    --p-maxhits all \
    --p-maxrejects all \
    --p-weak-id ${weak_id} \
    --p-maxaccepts ${maxaccepts} \
    --p-min-consensus 0.51 \
    --o-classification NERRS_18s_vsearch_taxonomy_10accepts_90perc-crux \
    --o-search-results NERRS_18s_vsearch_taxonomy_10accepts_90perc-tophits-crux
 
refreads=${refreads:-ref-dbs/SILVA/silva-138-99-seqs.qza}
reftax=${reftax:-ref-dbs/SILVA/silva-138-99-tax.qza}
sklearn=${sklearn:-ref-dbs/SILVA/silva-138-99-seqs-pid_0.65-classifier.qza}

qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/NERRS_18s_5_20_24_rep-seqs.qza \
  --i-classifier ${sklearn} \
  --i-reference-reads ${refreads} \
  --i-reference-taxonomy  ${reftax} \
  --p-threads 24 \
  --p-query-cov ${query_cov} \
  --p-perc-identity ${perc_identity} \
  --p-maxrejects all \
  --p-maxaccepts ${maxaccepts} \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva \
  --o-search-results NERRS_18s_vsearch_taxonomy_10accepts_90perc-tophits-silva

############################################################################################
############################################################################################

############################################################################################
############################################################################################
### Filter table
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_table.qza \
  --m-metadata-file results/qiime-swmp-sample-metadata.tsv \
  --o-filtered-table results/NERRS_18s_9_12_24_filtered-table.qza

qiime feature-table filter-features \
    --i-table results/NERRS_18s_5_20_24_filtered-table.qza \
    --m-metadata-file results/NERRS_18s_5_20_24_hybrid_taxonomy.qza \
    --o-filtered-table results/filtered-by-features_NERRS_18s_5_20_24.qza

### BARPLOT
qiime taxa barplot \
     --i-table results/filtered-by-features_NERRS_18s_5_20_24.qza \
     --i-taxonomy results/NERRS_18s_5_20_24_hybrid_taxonomy.qza \
     --o-visualization results/NERRS_18s_5_20_24_taxa-barplot \
     --m-metadata-file metadata/18s_sample_metadata_5-20.tsv

### Make output tables
conda activate qiime2R
export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
qiime2_output_tables.r results/NERRS_18s_5_20_24_filtered-table.qza results/NERRS_18s_5_20_24_hybrid_taxonomy.qza results/NERRS_18s_5_20_24_rep-seqs.qza NERRS_18s_5_20_24
############################################################################################
############################################################################################





############################################################################################
## BLANKS ##################################################################################
### Filter table
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_5_20_24_table.qza \
  --m-metadata-file metadata/18s_blank_metadata_5-20.tsv \
  --o-filtered-table results/NERRS_18s_5_20_24_blanks-table.qza

### Make output tables
conda activate qiime2R
export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
qiime2_output_tables.r results/NERRS_18s_5_20_24_blanks-table.qza results/NERRS_18s_5_20_24_hybrid_taxonomy.qza results/NERRS_18s_5_20_24_rep-seqs.qza NERRS_18s_5_20_24_blanks
############################################################################################
############################################################################################
