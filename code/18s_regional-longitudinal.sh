
## Gemelli tutorial
#qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
#qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
#qiime longitudinal feature-volatility


## Phylogenetic Robust Aitchison PCA (RP-PCA) with taxonomy
## qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
## https://github.com/biocore/gemelli/blob/8de57acf564f0abbf9accf9a9486a85404078414/ipynb/tutorials/Phylogenetic-RPCA-moving-pictures.ipynb#L113
for rgn in Gulf NE N-Pacific Pacific-Island; do 

  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table ${rgn}_NERRS_18s-table.qza \
      --i-phylogeny NERRS_18s_9_12_24_rooted-tree.qza \
      --m-taxonomy-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 10 \
      --p-min-sample-count 500 \
      --o-biplot regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-ordination.qza \
      --o-distance-matrix regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-distance.qza \
      --o-counts-by-node-tree regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-tree.qza \
      --o-counts-by-node regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-table.qza \
      --o-t2t-taxonomy regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-taxonomy.qza

  qiime empress community-plot\
    --i-tree regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-tree.qza\
    --i-feature-table regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-table.qza\
    --i-pcoa regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-ordination.qza\
    --m-sample-metadata-file qiime-swmp-corrected-sample-metadata.tsv\
    --m-feature-metadata-file regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-taxonomy.qza\
    --p-filter-missing-features\
    --p-number-of-features 50\
    --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_phylo-empress.qzv
done



## CTF with gemelli (longitudinal volatility (of PC1) interactive qiime plots)
## qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
## https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/IBD-Tutorial-QIIME2-CLI.md
## run the python custom script to get the feature volatility 
python /home/users/jtm1171/code/nerrs/18s_format-subject-metadata-gemelli.py qiime-swmp-corrected-sample-metadata.tsv

for rgn in Gulf NE N-Pacific Pacific-Island; do 
  ## Timepoint by site 
  qiime gemelli ctf\
      --i-table ${rgn}_NERRS_18s-table.qza\
      --m-sample-metadata-file qiime-swmp-corrected-sample-metadata.tsv \
      --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-state-column Quarter_num\
      --p-individual-id-column Site_Corrected\
      --output-dir regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf

  qiime longitudinal volatility \
      --m-metadata-file regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/state_subject_ordination.qza \
      --p-state-column Quarter_num \
      --p-individual-id-column subject_id \
      --p-default-group-column NERR \
      --p-default-metric PC1 \
      --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/${rgn}_state_subject_ordination.qzv

    qiime emperor biplot\
      --i-biplot regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/subject_biplot.qza \
      --m-sample-metadata-file  qiime-swmp-corrected-sample-metadata-subject-metadata.tsv \
      --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-number-of-features 20\
      --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/${rgn}_subject_biplot.qzv
done


## Longitudinal volitility
## https://docs.qiime2.org/2020.2/tutorials/longitudinal/

for rgn in Gulf NE N-Pacific Pacific-Island; do 

  qiime longitudinal feature-volatility \
    --i-table ${rgn}_NERRS_18s-table.qza \
    --m-metadata-file qiime-corrected-sample-metadata.tsv \
    --p-state-column Quarter_num \
    --p-individual-id-column Site_Corrected \
    --p-n-estimators 10 \
    --p-random-state 17 \
    --output-dir regional_core-diversity/${rgn}_with-repl/${rgn}_ecam-feat-volatility
done







