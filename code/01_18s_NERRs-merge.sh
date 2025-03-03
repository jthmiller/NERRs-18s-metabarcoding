## Merge 18s data

## NOTE: These samples ran twice: 
      2 AB161w0208232a
      2 AB181w0208231a
      2 GBGBw0213231
      2 GBLRw0213231
      2 HEKKw0202231
      2 HER9w0202231


## Rename samples specific to group. 
for run in runs/* ; do
echo $run 

    set=$( basename $run )
    tag=$( echo $set | sed 's/NERR//' | sed 's/-18SNX//' )
    echo -e "#SampleID\tSampleID-new" > results/${set}_rename.tsv
    ls $run/reads/poly-G-trimmed/ | cut -f1 -d'_' | sort | uniq | awk -v tag="$tag" -v OFS='\t' '{print $1,tag"_"$1}' >> results/${set}_rename.tsv

    qiime feature-table group \
      --i-table $run/qiime_out/${set}_table.qza \
      --p-axis sample \
      --m-metadata-file results/${set}_rename.tsv \
      --m-metadata-column SampleID-new \
      --p-mode sum \
      --o-grouped-table results/${set}_reindexed-table.qza
done 


for run in runs/NERRGBHEWE-18SNX091323 runs/NERRJCMAPBSF-18SNX091423 runs/NERRABGTM-18SNX091523 ; do
echo $run 
    set=$( basename $run )
    tag=$( echo $set | sed 's/NERR//' | sed 's/-18SNX//' )
    echo -e "#SampleID\tSampleID-new" > results/${set}_rename.tsv
    ls $run/reads/poly-G-trimmed/ | cut -f1 -d'_' | sort | uniq | awk -v tag="$tag" -v OFS='\t' '{print $1,tag"_"$1}' >> results/${set}_rename.tsv

    qiime feature-table group \
      --i-table $run/qiime_out/${set}_table.qza \
      --p-axis sample \
      --m-metadata-file results/${set}_rename.tsv \
      --m-metadata-column SampleID-new \
      --p-mode sum \
      --o-grouped-table results/${set}_reindexed-table.qza
done 

## ## added NERRSS-18SNX051723
## ## Merge vsearch taxonomy
## find runs -name *_vsearch_taxonomy.qza -exec echo "--i-data {} \\"  \; | grep -v tophit
## qiime feature-table merge-taxa \
##     --i-data runs/NERRSS-18SNX051723/qiime_out/NERRSS-18SNX051723_07252023_vsearch_taxonomy.qza \
##     --i-data runs/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_07032023_vsearch_taxonomy.qza \
##     --i-data runs/NERRGTM-18SNX052223/qiime_out/NERRGTM-18SNX052223_06272023_vsearch_taxonomy.qza \
##     --i-data runs/NERRJCWEMA-18SNX050923/qiime_out/NERRJCWEMA-18SNX050923_06272023_vsearch_taxonomy.qza \
##     --i-data runs/NERRSFPB-18SNX051723/qiime_out/NERRSFPB-18SNX051723_06272023_vsearch_taxonomy.qza \
##     --o-merged-data results/NERRS_18s_vsearch_taxonomy.qza
## 
## ## Merge vsearch taxonomy tophit
## 
## find runs -name *_vsearch_taxonomy.qza -exec echo "--i-data {} \\"  \; | grep tophit
## qiime feature-table merge-taxa \
##     --i-data runs/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_07242023_tophit_vsearch_taxonomy.qza \
##     --i-data runs/NERRGTM-18SNX052223/qiime_out/NERRGTM-18SNX052223_07242023_tophit_vsearch_taxonomy.qza \
##     --i-data runs/NERRJCWEMA-18SNX050923/qiime_out/NERRJCWEMA-18SNX050923_07242023_tophit_vsearch_taxonomy.qza \
##     --i-data runs/NERRSFPB-18SNX051723/qiime_out/NERRSFPB-18SNX051723_07242023_tophit_vsearch_taxonomy.qza \
##     --i-data runs/NERRSS-18SNX051723/qiime_out/NERRSS-18SNX051723_07242023_tophit_vsearch_taxonomy.qza \
##     --o-merged-data results/NERRS_18s_vsearch_tophit.qza
## 

## merge table
## find results -name *_reindexed-table.qza -exec echo "--i-tables {} \\"  \;
qiime feature-table merge \
    --i-tables results/NERRABGBHE-18SNX050923_reindexed-table.qza \
    --i-tables results/NERRGTM-18SNX052223_reindexed-table.qza \
    --i-tables results/NERRJCWEMA-18SNX050923_reindexed-table.qza \
    --i-tables results/NERRSFPB-18SNX051723_reindexed-table.qza \
    --i-tables results/NERRSS-18SNX051723_reindexed-table.qza \
    --i-tables results/NERRGBHEWE-18SNX091323_reindexed-table.qza \
    --i-tables results/NERRJCMAPBSF-18SNX091423_reindexed-table.qza \
    --i-tables results/NERRABGTM-18SNX091523_reindexed-table.qza \
    --o-merged-table results/NERRS_18s_reindexed-10-13_table.qza


## merge asvs
## find runs -name *rep-seqs.qza -exec echo "--i-data {} \\"  \;
qiime feature-table merge-seqs \
    --i-data runs/NERRABGBHE-18SNX050923/qiime_out/NERRABGBHE-18SNX050923_rep-seqs.qza \
    --i-data runs/NERRGTM-18SNX052223/qiime_out/NERRGTM-18SNX052223_rep-seqs.qza \
    --i-data runs/NERRJCWEMA-18SNX050923/qiime_out/NERRJCWEMA-18SNX050923_rep-seqs.qza \
    --i-data runs/NERRSFPB-18SNX051723/qiime_out/NERRSFPB-18SNX051723_rep-seqs.qza \
    --i-data runs/NERRSS-18SNX051723/qiime_out/NERRSS-18SNX051723_rep-seqs.qza \
    --i-data runs/NERRGBHEWE-18SNX091323/qiime_out/NERRGBHEWE-18SNX091323_rep-seqs.qza \
    --i-data runs/NERRJCMAPBSF-18SNX091423/qiime_out/NERRJCMAPBSF-18SNX091423_rep-seqs.qza \
    --i-data runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_rep-seqs.qza \
    --o-merged-data results/NERRS_18s_reindexed-10-13_rep-seqs.qza

head -n1 NERR-metadata_merge_samples.tsv > NEERs_18s_metadata.tsv
while read samp ; do 
eval "grep $( echo $samp | sed 's/18S//') NERR-metadata_merge_samples.tsv"
done <<< "$( cat ../results/*_rename.tsv | grep -v '^#' | cut -f1 )" | sort | uniq >> NEERs_18s_metadata.tsv

### metadata 
cp ../mifish/results/NERR-metadata-7_14.csv metadata/NERR-metadata-7_14.csv
tr ',' '\t' < metadata/NERR-metadata-7_14.csv > metadata/NERR-metadata-7_14.tsv
## cp metadata/NERR-metadata_merge_samples.tsv metadata/NERR-metadata_merge_samples.tsv



## regenerate metadata from mifish metadata
head -n1 metadata/NERR-metadata-7_14.tsv > metadata/NEERs_18s_10-13_metadata.tsv
while read samp ; do  
    oldname=$(echo $samp | cut -f1 -d' ' | sed 's/18S//' )
    newname=$(echo $samp | cut -f2 -d' ')
    grep $oldname metadata/NERR-metadata-7_14.tsv | awk -v OFS='\t' -v samp=$newname -v oldname=$oldname '$1==oldname{print samp,$0}'

done <<< "$( grep -v '^#' results/*_rename.tsv | cut -f2 -d':' )" >> metadata/NEERs_18s_10-13_metadata.tsv


qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_rep-seqs.qza \
   --o-alignment results/NERRS_18s_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_masked-aligned-rep-seqs.qza\
   --o-tree results/NERRS_18s_unrooted-tree.qza\
   --o-rooted-tree results/NERRS_18s_rooted-tree.qza\
   --p-n-threads 24


## Filter out samples not found
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_table.qza \
  --m-metadata-file metadata/NEERs_18s_metadata.tsv \
  --o-filtered-table results/NERRS_18s_filtered-table.qza

## Filter out controls
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_filtered-table.qza \
  --m-metadata-file metadata/NEERs_18s_metadata.tsv \
  --p-where "sample_blank='sample'" \
  --o-filtered-table results/NERRS_18s_samp_filtered-table.qza

#### Filter out AB
##qiime feature-table filter-samples \
##  --i-table results/NERRS_COI_samp_filtered-table.qza \
##  --m-metadata-file metadata/NEERs_COI_metadata.tsv \
##  --p-where "site='AB'" \
##  --p-exclude-ids \
##  --o-filtered-table results/NERRS_COI_NOAB_filtered-table.qza

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_rooted-tree.qza \
    --i-table results/NERRS_18s_samp_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/NEERs_18s_metadata.tsv \
    --output-dir results/all_sample_core-metrics-results

## ### Core Metrics
## qiime diversity core-metrics-phylogenetic \
##     --i-phylogeny results/NERRS_COI_rooted-tree.qza \
##     --i-table results/NERRS_COI_NOAB_filtered-table.qza \
##     --p-with-replacement \
##     --p-sampling-depth 500 \
##     --m-metadata-file metadata/NEERs_COI_metadata.tsv \
##     --output-dir results/no_AB_core-metrics-results

qiime feature-table summarize\
   --i-table results/NERRS_18s_filtered-table.qza \
   --m-sample-metadata-file metadata/NEERs_18s_metadata.tsv \
   --o-visualization results/NERRS_18s_filtered-table-summary

qiime taxa barplot \
     --i-table results/NERRS_18s_filtered-table.qza \
     --i-taxonomy results/NERRS_18s_vsearch_taxonomy.qza \
     --o-visualization results/NERRS_18s_filtered_taxa-barplot \
     --m-metadata-file metadata/NEERs_18s_metadata.tsv


# Classify rep seqs with pr2
qiime feature-classifier classify-sklearn \
    --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/pr2/pr2database-5.0.0/pr2_version_5.0.0_SSU_wPrimers-classifier.qza \
    --i-reads results/NERRS_18s_rep-seqs.qza \
    --o-classification results/.qza




### SILVA classifyer #### 
### Get non-protist taxa ids

# Classify rep seqs
qiime feature-classifier classify-sklearn \
--i-classifier /home/share/databases/SILVA_databases/silva-132-99-515-806-nb-classifier.qza \
--i-reads results/NERRS_18s_rep-seqs.qza \
--o-classification results/NERRS_18s_sklearn-taxonomy.qza

# Classify rep seqs
qiime feature-classifier classify-sklearn \
--i-classifier /home/share/databases/SILVA_databases/silva-132-99-515-806-nb-classifier.qza \
--i-reads results/NERRS_18s_reindexed-10-13_rep-seqs.qza \
--o-classification results/NERRS_18s_10-24_sklearn-taxonomy.qza
# Current primer
# fw=GTACACACCGCCCGTC	
# rv=TGATCCTTCTGCAGGTTCACCTAC
    
## old messed up primer
#fw=TGATCCTTCTGCAGGTTCACCTAC
#rv=GTACACACCGCCCGTC

qiime feature-classifier extract-reads \
  --i-sequences /home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-seqs.qza \
  --p-f-primer TGATCCTTCTGCAGGTTCACCTAC \
  --p-r-primer GTACACACCGCCCGTC \
  --p-n-jobs 12 \
  --o-reads 18s/silva/silva-138-99-seqs-extract-reads.qza 

qiime feature-classifier extract-reads \
  --i-sequences /home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-seqs.qza \
  --p-f-primer GTACACACCGCCCGTC \
  --p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
  --p-n-jobs 12 \
  --o-reads silva-138-99-seqs-extract-reads_former_primer_config.qza 



qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza   \
  --i-reference-taxonomy SILVA/silva-138-99-tax.qza \
  --o-classifier 18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 18s/silva/silva-138-99-seqs-extract-reads.qza   \
  --i-reference-taxonomy SILVA/silva-138-99-tax.qza \
  --o-classifier 18s/silva/silva-138-99_2022.8_nb-classifier.qza


ls -ltrh /home/unhAW/jtmiller/watts/ref-database/18s/silva/


cat \
  /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_export/dna-sequences.fasta \
  /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract_former_primer_config/dna-sequences.fasta \
  > /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract_former_primer_config/FW_RV_primers.fasta

qiime tools import \
  --input-path /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract_former_primer_config/FW_RV_primers.fasta \
  --output-path /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract_former_primer_config/FW_RV_18s_primer_seq.qza \
  --type 'FeatureData[Sequence]'


# Classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99_2022.8_nb-classifier.qza \
  --i-reads results/NERRS_18s_rep-seqs.qza \
  --o-classification results/NERRS_18s_sklearn-taxonomy.qza

# Classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99_2022.8_nb-classifier.qza \
  --i-reads results/NERRS_18s_10-13_rep-seqs.qza \
  --o-classification results/NERRS_18s_10-13_silva_sklearn-taxonomy.qza

# Classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99_2022.8_nb-classifier.qza \
  --i-reads results/NERRS_18s_reindexed-10-13_rep-seqs.qza \
  --o-classification results/NERRS_18s_10-24_sklearn-taxonomy.qza

# Classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza \
  --i-reads results/NERRS_18s_reindexed-10-13_rep-seqs.qza \
  --o-classification results/inverted_primer_NERRS_18s_10-24_sklearn-taxonomy.qza







nohup qiime feature-classifier classify-consensus-vsearch \
    --i-query results/NERRS_18s_reindexed-10-13_rep-seqs.qza \
    --i-reference-reads /home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza \
    --i-reference-taxonomy /home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-tax.qza \
    --p-maxaccepts all \
    --p-query-cov 0.8 \
    --p-perc-identity 0.90 \
    --o-classification results/silva-138-99-18s_vsearch_taxonomy \
    --p-threads 24 \
    --p-weak-id 0.8 & 


## # Classify rep seqs
## qiime feature-classifier classify-sklearn \
## --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/pr2/pr2database-5.0.0/pr2_version_5.0.0_SSU_wPrimers-classifier.qza \
## --i-reads results/NERRS_18s_rep-seqs.qza \
## --o-classification results/NERRS_18s_pr2_sklearn-taxonomy.qza

## # Classify rep seqs
## qiime feature-classifier classify-sklearn \
## --i-classifier /home/unhAW/jtmiller/watts/ref-database/18s/pr2/pr2database-5.0.0/pr2_version_5.0.0_SSU_wPrimers-classifier.qza \
## --i-reads results/NERRS_18s_10-13_rep-seqs.qza \
## --o-classification results/NERRS_18s_10-13_pr2_sklearn-taxonomy.qza

### Community weights?
# water-non-saline.qza
# water-saline.qza



## qiime feature-classifier extract-reads \
##  --i-sequences readytowear/data/silva_138_1/full_length/water-saline.qza \
##  --p-f-primer TGATCCTTCTGCAGGTTCACCTAC \
##  --p-r-primer GTACACACCGCCCGTC \
##  --p-n-jobs 12 \
##  --o-reads 18s/silva/water-saline-seqs-extract-reads.qza 
##
## qiime feature-classifier extract-reads \
##  --i-sequences readytowear/data/silva_138_1/full_length/water-saline.qza \
##  --p-f-primer TGATCCTTCTGCAGGTTCACCTAC \
##  --p-r-primer GTACACACCGCCCGTC \
##  --p-n-jobs 12 \
##  --o-reads 18s/silva/water-saline-seqs-extract-reads.qza 

qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_10-13_rep-seqs.qza \
   --o-alignment results/NERRS_18s_10-13_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_10-13_masked-aligned-rep-seqs.qza\
   --o-tree results/NERRS_18s_10-13_unrooted-tree.qza\
   --o-rooted-tree results/NERRS_18s_10-13_rooted-tree.qza\
   --p-n-threads 24


### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_10-13_rooted-tree.qza \
    --i-table results/NERRS_18s_10-13_table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/NEERs_18s_metadata.tsv \
    --output-dir results/NERRS_18s_10-13_core-metrics-results

#### Filter out samples not found
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_reindexed-10-13_table.qza \
  --m-metadata-file metadata/NEERs_18s_10-24_metadata.tsv \
  --o-filtered-table results/NEERs_18s_10-24_filtered-table.qza

## Filter out controls
qiime feature-table filter-samples \
  --i-table results/NEERs_18s_10-24_filtered-table.qza \
  --m-metadata-file metadata/NEERs_18s_10-24_metadata.tsv \
  --p-where "sample_blank='sample'" \
  --o-filtered-table results/NEERs_18s_10-24_samp_filtered-table.qza



#### Filter out AB
##qiime feature-table filter-samples \
##  --i-table results/NERRS_COI_samp_filtered-table.qza \
##  --m-metadata-file metadata/NEERs_COI_metadata.tsv \
##  --p-where "site='AB'" \
##  --p-exclude-ids \
##  --o-filtered-table results/NERRS_COI_NOAB_filtered-table.qza

qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_10-24_rep-seqs.qza \
   --o-alignment results/NERRS_18s_10-24_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_10-24_masked-aligned-rep-seqs.qza\
   --o-tree results/NERRS_18s_10-24_unrooted-tree.qza\
   --o-rooted-tree results/NERRS_18s_10-24_rooted-tree.qza\
   --p-n-threads 24

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_10-13_rooted-tree.qza \
    --i-table results/NEERs_18s_10-24_samp_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/NEERs_18s_10-24_metadata.tsv \
    --output-dir results/NEERs_18s_10-24_samp_filtered_core-metrics-results

## ### Core Metrics
## qiime diversity core-metrics-phylogenetic \
##     --i-phylogeny results/NERRS_COI_rooted-tree.qza \
##     --i-table results/NERRS_COI_NOAB_filtered-table.qza \
##     --p-with-replacement \
##     --p-sampling-depth 500 \
##     --m-metadata-file metadata/NEERs_COI_metadata.tsv \
##     --output-dir results/no_AB_core-metrics-results

qiime feature-table summarize\
   --i-table results/NEERs_18s_10-13_filtered-table.qza \
   --m-sample-metadata-file metadata/NEERs_18s_10-24_metadata.tsv \
   --o-visualization results/NEERs_18s_10-13_filtered-table-summary

qiime taxa barplot \
     --i-table results/NEERs_18s_10-13_samp_filtered-table.qza \
     --i-taxonomy results/NERRS_18s_10-13_silva_sklearn-taxonomy.qza \
     --o-visualization results/NEERs_18s_10-13_filtered_taxa-barplot \
     --m-metadata-file metadata/NEERs_18s_10-24_metadata.tsv

qiime2_output_tables.r qiime_out/${run}_table.qza qiime_out/${tag}_hybrid_taxonomy.qza qiime_out/${run}_rep-seqs.qza qiime_out/${tag}_ASV_table.csv qiime_out/${tag}_Species_table.csv









qiime feature-table merge \
    --i-tables runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_table.qza \
    --i-tables runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_table.qza \
    --i-tables runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_table.qza \
    --o-merged-table results/NERRS_18SNX012324_table.qza

## merge asvs
## find runs -name *rep-seqs.qza -exec echo "--i-data {} \\"  \;
qiime feature-table merge-seqs \
    --i-data runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_rep-seqs.qza \
    --o-merged-data results/NERRS_18SNX012324_rep-seqs.qza


refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza}
reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-tax.qza}
sklearn=${sklearn:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza}

    ## taxonomy
    maxaccepts=10
    query_cov=0.8 
    perc_identity=0.90 
    weak_id=0.80
    

qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/NERRS_18SNX012324_rep-seqs.qza \
  --i-classifier ${sklearn} \
  --i-reference-reads ${refreads} \
  --i-reference-taxonomy  ${reftax} \
  --p-threads 12 \
  --p-query-cov ${query_cov} \
  --p-perc-identity ${perc_identity} \
  --p-maxrejects all \
  --p-maxaccepts ${maxaccepts} \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/NERRS_18SNX012324_hybrid_taxonomy









#### Filter out samples not found
qiime feature-table filter-samples \
  --i-table results/NERRS_18SNX012324_table.qza \
  --m-metadata-file metadata/18s_NERRs_Nov2023.tsv \
  --o-filtered-table results/NERRS_18SNX012324_filtered-table.qza

## Filter out controls
qiime feature-table filter-samples \
  --i-table results/NERRS_18SNX012324_filtered-table.qza \
  --m-metadata-file metadata/18s_NERRs_Nov2023.tsv \
  --p-where "type='sample'" \
  --o-filtered-table results/NERRS_18SNX012324_samp_filtered-table.qza




qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18SNX012324_rep-seqs.qza \
   --o-alignment results/NERRS_18SNX012324_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18SNX012324_masked-aligned-rep-seqs.qza\
   --o-tree results/NERRS_18SNX012324_unrooted-tree.qza\
   --o-rooted-tree results/NERRS_18SNX012324_rooted-tree.qza\
   --p-n-threads 24

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18SNX012324_rooted-tree.qza \
    --i-table results/NERRS_18SNX012324_samp_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/18s_NERRs_Nov2023.tsv \
    --output-dir results/NERRS_18SNX012324_samp_filtered_core-metrics-results

qiime taxa barplot \
     --i-table results/NERRS_18SNX012324_samp_filtered-table.qza \
     --i-taxonomy results/NERRS_18SNX012324_hybrid_taxonomy.qza \
     --o-visualization results/NERRS_18SNX012324_filtered_taxa-barplot \
     --m-metadata-file metadata/18s_NERRs_Nov2023.tsv


qiime2_output_tables.r NERRS_18SNX012324_table.qza NERRS_18SNX012324_hybrid_taxonomy.qza NERRS_18SNX012324_rep-seqs.qza NERRS_18SNX012324_ASV_table.csv NERRS_18SNX012324_Species_table.csv



### Last two quarters
drwxr-xr-x 6 jtmiller unhaw   85 Mar  6 16:36 NERRGTMSSWE-18SNX012324
drwxr-xr-x 6 jtmiller unhaw   85 Mar  6 16:43 NERRHESF-18SNX012324
drwxr-xr-x 6 jtmiller unhaw   85 Mar  6 16:48 NERRMAPB-18SNX012324
drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:20 NERRGBGTMHE-18SNX013024
drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:38 NERRMAPBSFSS-18SNX013024
drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:38 NERRABJCWE-18SNX013024
NERRABGBJC-18SNX010424

### metadata 
tr ',' '\t' < metadata/18s_3_26_24.csv > metadata/18s_3_26_24.tsv


qiime feature-table merge-seqs \
    --i-data runs/NERRABGBJC-18SNX010424/qiime_out/NERRABGBJC-18SNX010424_rep-seqs.qza \
    --i-data runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_rep-seqs.qza \
    --i-data runs/NERRABJCWE-18SNX013024/qiime_out/NERRABJCWE-18SNX013024_rep-seqs.qza \
    --i-data runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRMAPBSFSS-18SNX013024/qiime_out/NERRMAPBSFSS-18SNX013024_rep-seqs.qza \
    --o-merged-data results/NERRS_18s_3_26_24_rep-seqs.qza


    ## taxonomy
    maxaccepts=10
    query_cov=0.8 
    perc_identity=0.90 
    weak_id=0.80
    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-tax.qza}
    sklearn=${sklearn:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza}


qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/NERRS_18s_3_26_24_rep-seqs.qza \
  --i-classifier ${sklearn} \
  --i-reference-reads ${refreads} \
  --i-reference-taxonomy  ${reftax} \
  --p-threads 12 \
  --p-query-cov ${query_cov} \
  --p-perc-identity ${perc_identity} \
  --p-maxrejects all \
  --p-maxaccepts ${maxaccepts} \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/NERRS_18s_3_26_24_hybrid_taxonomy
 
qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_3_26_24_rep-seqs.qza \
   --o-alignment results/NERRS_18s_3_26_24_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_3_26_24_masked-aligned-rep-seqs.qza \
   --o-tree results/NERRS_18s_3_26_24_unrooted-tree.qza \
   --o-rooted-tree results/NERRS_18s_3_26_24_rooted-tree.qza \
   --p-n-threads 24

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_3_26_24_rooted-tree.qza \
    --i-table results/NEERs_18s_03_26_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/18s_3_27_24.tsv \
    --output-dir results/NERRS_18s_3_26_24_samp_filtered_core-metrics-results




#### Filter out samples not found in lane that had extra samples in NERRMAPB-18SNX012324_table.qza
head -n1 metadata/18s_3_26_24.tsv > runs/NERRMAPB-18SNX012324/qiime_out/metadata
grep '^PB\|^MA' metadata/18s_3_26_24.tsv >> runs/NERRMAPB-18SNX012324/qiime_out/metadata
cut -f1-9,11-14 metadata/18s_3_26_24.tsv > metadata/18s_3_27_24.tsv

qiime feature-table filter-samples \
  --i-table runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_table.qza \
  --m-metadata-file runs/NERRMAPB-18SNX012324/qiime_out/metadata2 \
  --o-filtered-table runs/NERRMAPB-18SNX012324/qiime_out/filt_NERRMAPB-18SNX012324-table.qza

qiime feature-table merge \
    --i-tables runs/NERRABGBJC-18SNX010424/qiime_out/NERRABGBJC-18SNX010424_table.qza \
    --i-tables runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_table.qza \
    --i-tables runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_table.qza \
    --i-tables runs/NERRABJCWE-18SNX013024/qiime_out/NERRABJCWE-18SNX013024_table.qza \
    --i-tables runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_table.qza \
    ## --i-tables runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_table.qza \
    --i-tables runs/NERRMAPBSFSS-18SNX013024/qiime_out/NERRMAPBSFSS-18SNX013024_table.qza \
    --i-tables runs/NERRMAPB-18SNX012324/qiime_out/filt_NERRMAPB-18SNX012324-table.qza \
    --o-merged-table results/NERRS_18s_3_26_24_table.qza

#### Filter out samples not found
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_3_26_24_table.qza \
  --m-metadata-file metadata/18s_3_27_24.tsv \
  --o-filtered-table results/NEERs_18s_03_26_filtered-table.qza

## ## Filter out controls
## qiime feature-table filter-samples \
##   --i-table results/NEERs_18s_03_26_filtered-table.qza \
##   --m-metadata-file metadata/18s_3_26_24.tsv \
##   --p-where "type='sample'" \
##   --o-filtered-table results/NEERs_18s_03_26_samp_filtered-table.qza

qiime feature-table filter-features \
    --i-table results/NEERs_18s_03_26_filtered-table.qza \
    --m-metadata-file metadata/18s_3_27_24.tsv \
    --o-filtered-table results/NEERs_18s_03_26_ASV_filtered-table

qiime taxa barplot \
     --i-table  results/NEERs_18s_03_26_ASV_filtered-table.qza \
     --i-taxonomy results/NERRS_18s_3_26_24_hybrid_taxonomy.qza \
     --o-visualization results/NEERs_18s_03_26_filtered_taxa-barplot \
     --m-metadata-file metadata/18s_3_27_24.tsv

qiime2_output_tables.r \
  results/NEERs_18s_03_26_filtered-table.qza \
  results/NERRS_18s_3_26_24_hybrid_taxonomy.qza \
  results/NERRS_18s_3_26_24_rep-seqs.qza \
  results/NERRS_18s_3_26_24_ASV_table.csv \
  results/NERRS_18s_3_26_24_Species_table.csv




qiime diversity alpha-phylogenetic \
  --i-table results/NEERs_18s_03_26_filtered-table.qza \
  --i-phylogeny results/NERRS_18s_3_26_24_rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity results/NERRS_18s_3_26_24_samp_filtered_core-metrics-results/faith_pd


qiime diversity alpha-rarefaction \
    --i-table results/NEERs_18s_03_26_filtered-table.qza \
    --i-phylogeny results/NERRS_18s_3_26_24_rooted-tree.qza \
    --p-max-depth 300000 \
    --m-metadata-file metadata/18s_3_27_24.tsv \
    --p-min-depth 500 \
    --p-steps 15 \
    --o-visualization results/NERRS_18s_3_26_24_samp_filtered_core-metrics-results/alpha-rarefaction

qiime diversity alpha-group-significance \
    --i-alpha-diversity results/NERRS_18s_3_26_24_samp_filtered_core-metrics-results/faith_pd.qza \
    --m-metadata-file metadata/18s_3_27_24.tsv \
    --o-visualization results/NERRS_18s_3_26_24_samp_filtered_core-metrics-results/alpha-group-significance





ut are not present in the reference taxonomy: KY399207.1, FJ383531.1, FJ383805.1, MH449202.1

grep 'KY399207.1' klymus.fasta-reads-retained-tax.txt
grep 'KY399207.1' klymus.fasta-reads-retained2/

qiime feature-table filter-features \
    --i-table klymus.fasta-reads.qza \
    --m-metadata-file klymus.fasta-reads-retained-tax.qza \
    --o-filtered-table filter

398579 in klymus.fasta-reads-retained










##### 18s merge 5-9-24

qiime feature-table merge-seqs \
    --i-data runs/NERRGBGTMSFSS-18SNX032624/qiime_out/NERRGBGTMSFSS-18SNX032624_rep-seqs.qza \
    --i-data runs/NERRSFWE-18SNX030524/qiime_out/NERRSFWE-18SNX030524_rep-seqs.qza \
    --i-data runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_rep-seqs.qza \
    --i-data runs/NERRABGBHE-18SNX030524/qiime_out/NERRABGBHE-18SNX030524_rep-seqs.qza \
    --i-data runs/NERRABPB-18SNX022224/qiime_out/NERRABPB-18SNX022224_rep-seqs.qza \
    --i-data runs/NERRHEJCMAWE-18SNX031424/qiime_out/NERRHEJCMAWE-18SNX031424_rep-seqs.qza \
    --i-data runs/NERRJCMAPB-18SNX030524/qiime_out/NERRJCMAPB-18SNX030524_rep-seqs.qza \
    --o-merged-data results/NERRS_18s_5_9_24_rep-seqs.qza


qiime feature-table merge \
    --i-tables runs/NERRGBGTMSFSS-18SNX032624/qiime_out/NERRGBGTMSFSS-18SNX032624_table.qza \
    --i-tables runs/NERRSFWE-18SNX030524/qiime_out/NERRSFWE-18SNX030524_table.qza \
    --i-tables runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_table.qza \
    --i-tables runs/NERRABGBHE-18SNX030524/qiime_out/NERRABGBHE-18SNX030524_table.qza \
    --i-tables runs/NERRABPB-18SNX022224/qiime_out/NERRABPB-18SNX022224_table.qza \
    --i-tables runs/NERRHEJCMAWE-18SNX031424/qiime_out/NERRHEJCMAWE-18SNX031424_table.qza \
    --i-tables runs/NERRJCMAPB-18SNX030524/qiime_out/NERRJCMAPB-18SNX030524_table.qza \
    --o-merged-table results/NERRS_18s_5_9_24_table.qza


     43 NERRABGBJC-18SNX010424
     12 NERRABGTM-18SNX091523
     44 NERRABJCWE-18SNX013024
     42 NERRGTMSSWE-18SNX012324
     31 NERRHESF-18SNX012324
     23 NERRMAPB-18SNX012324
     57 NERRMAPBSFSS-18SNX013024






missing:
NERRABJCWE-18SNX013024



drwxr-xr-x 6 jtmiller unhaw   85 Jul  3  2023 NERRABGBHE-18SNX050923
drwxr-xr-x 6 jtmiller unhaw   85 Oct  4  2023 NERRGBHEWE-18SNX091323
drwxr-xr-x 6 jtmiller unhaw   85 Mar 26 15:39 NERRABGBJC-18SNX010424

while read sample ; do 

find . -maxdepth 3 -name *${sample}*R1_001.fastq.gz | cut -f2 -d'/' | sort | uniq 

done < missing.samples | sort | uniq -c






drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:20 NERRGBGTMHE-18SNX013024
drwxr-xr-x 6 jtmiller unhaw   85 May  9 09:56 NERRGBGTMSFSS-18SNX032624
drwxr-xr-x 6 jtmiller unhaw   85 May  9 10:41 NERRABGBHE-18SNX030524



?? 
drwxr-xr-x 6 jtmiller unhaw   85 Mar  6 16:36 NERRGTMSSWE-18SNX012324
drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:38 NERRMAPBSFSS-18SNX013024

drwxr-xr-x 6 jtmiller unhaw   85 Jun 26  2023 NERRSS-18SNX051723
drwxr-xr-x 6 jtmiller unhaw   85 Mar  6 16:36 NERRGTMSSWE-18SNX012324
drwxr-xr-x 6 jtmiller unhaw   85 Mar 15 11:38 NERRMAPBSFSS-18SNX013024






    ## taxonomy
    maxaccepts=10
    query_cov=0.8 
    perc_identity=0.90 
    weak_id=0.80
    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-tax.qza}
    sklearn=${sklearn:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza}


qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query results/NERRS_18s_5_9_24_rep-seqs.qza \
  --i-classifier ${sklearn} \
  --i-reference-reads ${refreads} \
  --i-reference-taxonomy  ${reftax} \
  --p-threads 12 \
  --p-query-cov ${query_cov} \
  --p-perc-identity ${perc_identity} \
  --p-maxrejects all \
  --p-maxaccepts ${maxaccepts} \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification results/NERRS_18s_5_9_24_hybrid_taxonomy
 
qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_5_9_24_rep-seqs.qza \
   --o-alignment results/NERRS_18s_5_9_24_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_5_9_24_masked-aligned-rep-seqs.qza \
   --o-tree results/NERRS_18s_5_9_24_unrooted-tree.qza \
   --o-rooted-tree results/NERRS_18s_5_9_24_rooted-tree.qza \
   --p-n-threads 24

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_5_9_24_rooted-tree.qza \
    --i-table results/NEERs_18s_03_26_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/18s_3_27_24.tsv \
    --output-dir results/NERRS_18s_5_9_24_samp_filtered_core-metrics-results


qiime tools export --input-path results/NERRS_18s_5_9_24_table.qza
biom convert -i results/NERRS_18s_5_9_24_table/feature-table.biom -o results/NERRS_18s_5_9_24_table/feature-table.tsv --to-tsv




results/NERRS_18s_20_9_24_table.qza

############ REMERGE


while read sample ; do 

find . -maxdepth 3 -name *${sample}*R1_001.fastq.gz | cut -f2 -d'/' | sort | uniq 
echo  ' '

done < missing.samples | sort | uniq -c

## found in multip[le runs 

for samp in AB221w0523232 AB161w0523232 AB181tb052323 AB151w0523233 AB181w0523233 AB221w0523233 AB221w0523231 AB181w0523232 ; do
find . -maxdepth 3 -name *${samp}*R1_001.fastq.gz | cut -f2 -d'/' | sort | uniq
done | sort | uniq -c

NERRABGTM-18SNX091523
NERRABGBHE-18SNX030524


     NERRABGTM-18SNX091523

     24 NERRABGTM-18SNX091523
     58 NERRGBHEWE-18SNX091323
     79 NERRJCMAPBSF-18SNX091423

     55 NERRJCMAPB-18SNX030524
     66 NERRABGBHE-18SNX030524
     40 NERRSFWE-18SNX030524


NERRABGBHE-18SNX030524


runs/NERRABGBHE-18SNX030524 (this is the rerun)
runs/NERRABGTM-18SNX091523 (of these samples)

##### 18s merge 5-20-24

qiime feature-table merge-seqs \
    --i-data runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_rep-seqs.qza \
    --i-data runs/NERRGBGTMSFSS-18SNX032624/qiime_out/NERRGBGTMSFSS-18SNX032624_rep-seqs.qza \
    --i-data runs/NERRSFWE-18SNX030524/qiime_out/NERRSFWE-18SNX030524_rep-seqs.qza \
    --i-data runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_rep-seqs.qza \
    --i-data runs/NERRABGBHE-18SNX030524/qiime_out/NERRABGBHE-18SNX030524_rep-seqs.qza \
    --i-data runs/NERRABPB-18SNX022224/qiime_out/NERRABPB-18SNX022224_rep-seqs.qza \
    --i-data runs/NERRHEJCMAWE-18SNX031424/qiime_out/NERRHEJCMAWE-18SNX031424_rep-seqs.qza \
    --i-data runs/NERRJCMAPB-18SNX030524/qiime_out/NERRJCMAPB-18SNX030524_rep-seqs.qza \
    --i-data runs/NERRABGBJC-18SNX010424/qiime_out/NERRABGBJC-18SNX010424_rep-seqs.qza \
    --i-data runs/NERRABJCWE-18SNX013024/qiime_out/NERRABJCWE-18SNX013024_rep-seqs.qza \
    --i-data runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_rep-seqs.qza \
    --i-data runs/NERRMAPBSFSS-18SNX013024/qiime_out/NERRMAPBSFSS-18SNX013024_rep-seqs.qza \
    --o-merged-data  ~/NERRs-results/18s/NERRS_18s_5_20_24_rep-seqs.qza

#   --i-tables runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_table.qza \
qiime feature-table merge \
    --i-tables ~/NERRs-results/18s/NERRABGTM-18SNX091523_filtered-table.qza \
    --i-tables runs/NERRGBGTMSFSS-18SNX032624/qiime_out/NERRGBGTMSFSS-18SNX032624_table.qza \
    --i-tables runs/NERRSFWE-18SNX030524/qiime_out/NERRSFWE-18SNX030524_table.qza \
    --i-tables runs/NERRGBGTMHE-18SNX013024/qiime_out/NERRGBGTMHE-18SNX013024_table.qza \
    --i-tables runs/NERRABGBHE-18SNX030524/qiime_out/NERRABGBHE-18SNX030524_table.qza \
    --i-tables runs/NERRABPB-18SNX022224/qiime_out/NERRABPB-18SNX022224_table.qza \
    --i-tables runs/NERRHEJCMAWE-18SNX031424/qiime_out/NERRHEJCMAWE-18SNX031424_table.qza \
    --i-tables runs/NERRJCMAPB-18SNX030524/qiime_out/NERRJCMAPB-18SNX030524_table.qza \
    --i-tables runs/NERRABGBJC-18SNX010424/qiime_out/NERRABGBJC-18SNX010424_table.qza \
    --i-tables runs/NERRABJCWE-18SNX013024/qiime_out/NERRABJCWE-18SNX013024_table.qza \
    --i-tables runs/NERRGTMSSWE-18SNX012324/qiime_out/NERRGTMSSWE-18SNX012324_table.qza \
    --i-tables runs/NERRHESF-18SNX012324/qiime_out/NERRHESF-18SNX012324_table.qza \
    --i-tables runs/NERRMAPB-18SNX012324/qiime_out/NERRMAPB-18SNX012324_table.qza \
    --i-tables runs/NERRMAPBSFSS-18SNX013024/qiime_out/NERRMAPBSFSS-18SNX013024_table.qza \
    --o-merged-table  ~/NERRs-results/18s/NERRS_18s_table.qza


qiime tools export --input-path results/NERRS_18s_20_9_24_table.qza --output-path results/NERRS_18s_20_9_24_table
biom convert -i results/NERRS_18s_20_9_24_table/feature-table.biom -o results/NERRS_18s_20_9_24_table/feature-table.tsv --to-tsv
head -n2 results/NERRS_18s_20_9_24_table/feature-table.tsv > samps


# runs.txt, add to metadata
NERRGBGTMSFSS-18SNX032624
NERRSFWE-18SNX030524
NERRGBGTMHE-18SNX013024
NERRABGBHE-18SNX030524
NERRABPB-18SNX022224
NERRHEJCMAWE-18SNX031424
NERRJCMAPB-18SNX030524
NERRABGBJC-18SNX010424
NERRABJCWE-18SNX013024
NERRGTMSSWE-18SNX012324
NERRHESF-18SNX012324
NERRMAPB-18SNX012324
NERRMAPBSFSS-18SNX013024

NERRABGTM-18SNX091523

GTM samples did not get a re-run when AB did?

while read run ; do
    qiime metadata tabulate \
        --m-input-file runs/${run}/qiime_out/${run}_dns.qza \
        --o-visualization runs/${run}/qiime_out/${run}_dns
    
    qiime tools export \
        --input-path runs/${run}/qiime_out/${run}_dns.qzv \
        --output-path runs/${run}/qiime_out/${run}_dns_export
done < runs.txt






## taxonomy
maxaccepts=10
query_cov=0.8 
perc_identity=0.90 
weak_id=0.80
refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/silva-138-99-seqs-extract-reads_former_primer_config.qza}
reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99-tax.qza}
sklearn=${sklearn:-/home/unhAW/jtmiller/watts/ref-database/18s/silva/inverted_primer_silva-138-99_2022.8_nb-classifier.qza}

qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query ~/NERRs-results/18s/NERRS_18s_5_20_24_rep-seqs.qza \
  --i-classifier ${sklearn} \
  --i-reference-reads ${refreads} \
  --i-reference-taxonomy  ${reftax} \
  --p-threads 12 \
  --p-query-cov ${query_cov} \
  --p-perc-identity ${perc_identity} \
  --p-maxrejects all \
  --p-maxaccepts ${maxaccepts} \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0.7 \
  --o-classification ~/NERRs-results/18s/NERRS_18s_5_20_24_hybrid_taxonomy
 
qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences results/NERRS_18s_5_20_24_rep-seqs.qza \
   --o-alignment results/NERRS_18s_5_20_24_aligned-rep-seqs \
   --o-masked-alignment results/NERRS_18s_5_20_24_masked-aligned-rep-seqs.qza \
   --o-tree results/NERRS_18s_5_20_24_unrooted-tree.qza \
   --o-rooted-tree results/NERRS_18s_5_20_24_rooted-tree.qza \
   --p-n-threads 24

### Filter table
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_5_20_24_table.qza \
  --m-metadata-file metadata/18s_sample_metadata_5-20.tsv \
  --o-filtered-table results/NERRS_18s_5_20_24_filtered-table.qza

### Core Metrics
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
    --i-table results/NERRS_18s_5_20_24_filtered-table.qza \
    --p-with-replacement \
    --p-sampling-depth 500 \
    --m-metadata-file metadata/18s_sample_metadata_5-20.tsv \
    --output-dir results/NERRS_18s_5_20_24_samp_filtered_core-metrics-results

qiime diversity alpha-phylogenetic \
  --i-table results/NERRS_18s_5_20_24_filtered-table.qza \
  --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity results/NERRS_18s_5_20_24_samp_filtered_core-metrics-results/faith_pd

qiime diversity alpha-rarefaction \
    --i-table results/NERRS_18s_5_20_24_filtered-table.qza \
    --i-phylogeny results/NERRS_18s_5_20_24_rooted-tree.qza \
    --p-max-depth 300000 \
    --m-metadata-file metadata/18s_sample_metadata_5-20.tsv \
    --p-min-depth 1000 \
    --p-steps 15 \
    --o-visualization results/NERRS_18s_5_20_24_samp_filtered_core-metrics-results/alpha-rarefaction

qiime diversity alpha-group-significance \
    --i-alpha-diversity results/NERRS_18s_5_20_24_samp_filtered_core-metrics-results/faith_pd.qza \
    --m-metadata-file metadata/18s_sample_metadata_5-20.tsv \
    --o-visualization results/NERRS_18s_5_20_24_samp_filtered_core-metrics-results/alpha-group-significance


qiime feature-table filter-features \
    --i-table results/NERRS_18s_5_20_24_filtered-table.qza \
    --m-metadata-file results/NERRS_18s_5_20_24_hybrid_taxonomy.qza \
    --o-filtered-table results/filtered-by-features_NERRS_18s_5_20_24.qza

qiime taxa barplot \
     --i-table results/filtered-by-features_NERRS_18s_5_20_24.qza \
     --i-taxonomy results/NERRS_18s_5_20_24_hybrid_taxonomy.qza \
     --o-visualization results/NERRS_18s_5_20_24_taxa-barplot \
     --m-metadata-file metadata/18s_sample_metadata_5-20.tsv

### Make output tables
conda activate qiime2R
export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
qiime2_output_tables.r results/NERRS_18s_5_20_24_filtered-table.qza results/NERRS_18s_5_20_24_hybrid_taxonomy.qza results/NERRS_18s_5_20_24_rep-seqs.qza NERRS_18s_5_20_24




## BLANKS
### Filter table
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_5_20_24_table.qza \
  --m-metadata-file metadata/18s_blank_metadata_5-20.tsv \
  --o-filtered-table results/NERRS_18s_5_20_24_blanks-table.qza

### Make output tables
conda activate qiime2R
export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
qiime2_output_tables.r results/NERRS_18s_5_20_24_blanks-table.qza results/NERRS_18s_5_20_24_hybrid_taxonomy.qza results/NERRS_18s_5_20_24_rep-seqs.qza NERRS_18s_5_20_24_blanks







## Try CRUX database
## taxonomy
maxaccepts=10
query_cov=0.8 
perc_identity=0.90 
weak_id=0.80
refreads=${refreads:-/home/users/jtm1171/old-home/watts/ref-database/18s/CRUX/18S_.fasta.qza}
reftax=${reftax:-/home/users/jtm1171/old-home/watts/ref-database/18s/CRUX/18S_taxonomy.qza}
sklearn=${sklearn:-/home/users/jtm1171/old-home/watts/ref-database/18s/CRUX/18s_crux_nb-classifier.qza}

qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query ~/NERRs-results/18s/NERRS_18s_5_20_24_rep-seqs.qza \
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
  --o-classification ~/NERRs-results/18s/NERRS_18s_CRUX_hybrid_taxonomy
 


 ### Filter table
qiime feature-table filter-samples \
  --i-table /home-wd/home/unhAW/jtmiller/watts/data/NERR/18s/runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_table.qza \
  --m-metadata-file ~/NERRs-results/18s/filter.GTM.metadata \
  --o-filtered-table ~/NERRs-results/18s/NERRABGTM-18SNX091523_filtered-table.qza




  ### Filter table
qiime feature-table filter-samples \
  --i-table ~/NERRs-results/18s/NERRS_18s_table.qza \
  --m-metadata-file ~/NERRs-results/18s/qiime-swmp-sample-metadata.tsv \
  --o-filtered-table ~/NERRs-results/18s/NERRS_18s_9_12_24_filtered-table.qza


qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences NERRS_18s_5_20_24_rep-seqs.qza \
   --o-alignment NERRS_18s_5_20_24_aligned-rep-seqs \
   --o-masked-alignment NERRS_18s_5_20_24_masked-aligned-rep-seqs.qza \
   --o-tree NERRS_18s_5_20_24_unrooted-tree.qza \
   --o-rooted-tree NERRS_18s_5_20_24_rooted-tree.qza \
   --p-n-threads 24





## Try CRUX database
## taxonomy
maxaccepts=10
query_cov=0.8 
perc_identity=0.90 
weak_id=0.80
refreads=${refreads:-/home/users/jtm1171/refdbs/18s/CRUX/18S_.fasta.qza}
reftax=${reftax:-/home/users/jtm1171/refdbs/18s/CRUX/18S_taxonomy.qza}
sklearn=${sklearn:-/home/users/jtm1171/refdbs/18s/CRUX/18s_crux_nb-classifier.qza}

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
    --o-classification NERRS_18s_vsearch_taxonomy_10accepts_90perc \
    --o-search-results NERRS_18s_vsearch_taxonomy_10accepts_90perc-tophits

 


## Try CRUX database
## taxonomy
maxaccepts=10
query_cov=0.8 
perc_identity=0.90 
weak_id=0.80
refreads=${refreads:-/home/users/jtm1171/refdbs/18s/SILVA/silva-138-99-seqs.qza}
reftax=${reftax:-/home/users/jtm1171/refdbs/18s/SILVA/silva-138-99-tax.qza}
sklearn=${sklearn:-/home/users/jtm1171/refdbs/18s/SILVA/silva-138-99-seqs-pid_0.65-classifier.qza}

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
    --o-classification NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva \
    --o-search-results NERRS_18s_vsearch_taxonomy_10accepts_90perc-tophits-silva
