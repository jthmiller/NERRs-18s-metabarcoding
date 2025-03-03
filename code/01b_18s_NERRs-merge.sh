
##### 18s merge 5-20-24
 ### Filter table
qiime feature-table filter-samples \
  --i-table runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_table.qza \
  --m-metadata-file runs/NERRABGTM-18SNX091523/qiime_out/filter.GTM.metadata \
  --o-filtered-table runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_filtered-table.qza

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
    --o-merged-data ~/NERRs-results/18s/NERRS_18s_5_20_24_rep-seqs.qza

#   --i-tables runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_table.qza \
qiime feature-table merge \
    --i-tables runs/NERRABGTM-18SNX091523/qiime_out/NERRABGTM-18SNX091523_filtered-table.qza \
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
    --o-merged-table results/NERRS_18s_table.qza


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

### Tables for adding sample read counts to metadata
while read run ; do
    qiime metadata tabulate \
        --m-input-file runs/${run}/qiime_out/${run}_dns.qza \
        --o-visualization runs/${run}/qiime_out/${run}_dns
    
    qiime tools export \
        --input-path runs/${run}/qiime_out/${run}_dns.qzv \
        --output-path runs/${run}/qiime_out/${run}_dns_export
done < runs.txt




for file in $( ls ) ; do
    echo $file
done 



while read file ; do 
    echo $file
done < listfiles.txt  


