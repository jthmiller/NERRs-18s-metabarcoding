


for rgn in Gulf NE N-Pacific Pacific-Island; do 
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
        --i-table ${rgn}_NERRS_18s-table.qza \
        --p-with-replacement \
        --p-sampling-depth 1000 \
        --m-metadata-file qiime-corrected-sample-metadata-region.tsv \
        --output-dir regional_core-diversity/${rgn}_with-repl/core-metrics-phylogenetic

    ## Alpha diversity
    qiime diversity alpha-phylogenetic \
      --i-table ${rgn}_NERRS_18s-table.qza \
      --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
      --p-metric faith_pd \
      --o-alpha-diversity regional_core-diversity/${rgn}_with-repl/${rgn}_faith_pd

    qiime diversity alpha-rarefaction \
        --i-table ${rgn}_NERRS_18s-table.qza \
        --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
        --p-max-depth 150000 \
        --m-metadata-file qiime-corrected-sample-metadata-region.tsv  \
        --p-min-depth 100 \
        --p-steps 15 \
        --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_alpha-rarefaction

    qiime diversity alpha-group-significance \
        --i-alpha-diversity regional_core-diversity/${rgn}_with-repl/${rgn}_faith_pd.qza \
        --m-metadata-file qiime-corrected-sample-metadata-region.tsv  \
        --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_alpha-group-significance
done 



### Beta diversity
for rgn in Gulf NE N-Pacific Pacific-Island; do 
    qiime diversity beta-phylogenetic \
        --i-table ${rgn}_NERRS_18s-table.qza \
        --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
        --p-metric unweighted_unifrac \
        --o-distance-matrix regional_core-diversity/${rgn}_with-repl/${rgn}_unweighted_unifrac_distance_matrix.qza

    qiime diversity beta-phylogenetic \
        --i-table ${rgn}_NERRS_18s-table.qza \
        --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
        --p-metric weighted_unifrac \
        --o-distance-matrix regional_core-diversity/${rgn}_with-repl/${rgn}_weighted_unifrac_distance_matrix.qza

  for parm in salinity Quarter_txt NERR ; do

      qiime diversity beta-group-significance \
          --i-distance-matrix regional_core-diversity/${rgn}_with-repl/${rgn}_unweighted_unifrac_distance_matrix.qza \
          --m-metadata-file qiime-corrected-sample-metadata-region.tsv  \
          --m-metadata-column ${parm} \
          --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_${parm}_beta-group-significance

      qiime diversity beta-group-significance \
          --i-distance-matrix regional_core-diversity/${rgn}_with-repl/${rgn}_weighted_unifrac_distance_matrix.qza \
          --m-metadata-file qiime-corrected-sample-metadata-region.tsv  \
          --m-metadata-column ${parm} \
          --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_${parm}_weighted_beta-group-significance
  done

done

