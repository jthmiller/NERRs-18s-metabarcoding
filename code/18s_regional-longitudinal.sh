
## Gemelli tutorial
#qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
#qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
#qiime longitudinal feature-volatility


############################################################################################
############################################################################################
## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "NERR='HE'" \
  --o-filtered-table results/Pacific-Island_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='SE'" \
  --o-filtered-table results/SE_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='NE'" \
  --o-filtered-table results/NE_NERRS_18s-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='N-Pacific'" \
  --o-filtered-table results/N-Pacific_NERRS_18s-table.qza
############################################################################################
############################################################################################

############################################################################################
############################################################################################
## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "NERR='HE'" \
  --o-filtered-table results/Pacific-Island_genus-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='SE'" \
  --o-filtered-table results/SE_genus-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='NE'" \
  --o-filtered-table results/NE_genus-table.qza

## Filter by region
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_genus-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region='N-Pacific'" \
  --o-filtered-table results/N-Pacific_genus-table.qza
############################################################################################
############################################################################################

############################################################################################
############################################################################################
## Phylogenetic Robust Aitchison PCA (RP-PCA) with taxonomy
## qiime gemelli phylogenetic-rpca-with-taxonomy -> qiime empress community-plot
## https://github.com/biocore/gemelli/blob/8de57acf564f0abbf9accf9a9486a85404078414/ipynb/tutorials/Phylogenetic-RPCA-moving-pictures.ipynb#L113
for rgn in SE NE N-Pacific Pacific-Island; do 

  qiime gemelli phylogenetic-rpca-with-taxonomy \
      --i-table results/${rgn}_NERRS_18s-table.qza \
      --i-phylogeny results/NERRS_18s_9_12_24_rooted-tree.qza \
      --m-taxonomy-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-min-feature-count 10 \
      --p-min-sample-count 500 \
      --p-min-feature-frequency 0.05 \
      --o-biplot results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-ordination.qza \
      --o-distance-matrix results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-distance.qza \
      --o-counts-by-node-tree results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-tree.qza \
      --o-counts-by-node results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-table.qza \
      --o-t2t-taxonomy results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-taxonomy.qza

  qiime empress community-plot\
    --i-tree results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-tree.qza\
    --i-feature-table results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-table.qza\
    --i-pcoa results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-ordination.qza\
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-taxonomy.qza\
    --p-filter-missing-features\
    --p-number-of-features 50\
    --o-visualization results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-empress.qzv

qiime qurro loading-plot \
    --i-ranks results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-ordination.qza \
    --i-table results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-taxonomy.qza \
    --o-visualization results/regional-phylogenetic-rpca-with-taxonomy/${rgn}_phylo-qurro_plot.qzv
done

############################################################################################
############################################################################################
cp results/regional-phylogenetic-rpca-with-taxonomy/*qzv ../jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/


for rgn in SE NE N-Pacific Pacific-Island; do
echo -e "
[${rgn}_phylo-empress](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/${rgn}_phylo-empress.qzv)
[${rgn}_qurro-phylogenetic-rpca-with-taxonomy](https://view.qiime2.org/visualization/?src=https://jthmiller.github.io/files/results/nerrs/core-metrics-results/phylogenetic/${rgn}_phylo-qurro_plot.qzv)
"
done


## CTF with gemelli (longitudinal volatility (of PC1) interactive qiime plots)
## qiime gemelli ctf -> qiime longitudinal volatility -> qiime emperor biplot
## https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/IBD-Tutorial-QIIME2-CLI.md
## run the python custom script to get the feature volatility 
python /home/users/jtm1171/code/nerrs/18s_format-subject-metadata-gemelli.py qiime-swmp-corrected-sample-metadata.tsv

for rgn in SE NE N-Pacific Pacific-Island; do 
  ## Timepoint by site 
  qiime gemelli ctf\
    --i-table results/${rgn}_NERRS_18s-table.qza\
    --m-sample-metadata-file metadata.tsv \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --p-min-feature-frequency 0.05 \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --o-subject-biplot results/regional-ctf/${rgn}-ctf_subject_biplot.qza \
    --o-state-biplot results/regional-ctf/${rgn}-ctf_state_biplot.qza \
    --o-distance-matrix results/regional-ctf/${rgn}-ctf_distance_matrix.qza \
    --o-state-subject-ordination results/regional-ctf/${rgn}-ctf_state-subject-ordination.qza \
    --o-state-feature-ordination results/regional-ctf/${rgn}-ctf_state-feature-ordination.qza \

qiime qurro loading-plot \
    --i-ranks results/regional-ctf/${rgn}-ctf_subject_biplot.qza \
    --i-table results/${rgn}_NERRS_18s-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --o-visualization results/regional-ctf/${rgn}_ctf-qurro_plot.qzv
done 



for rgn in SE NE N-Pacific Pacific-Island; do 

  qiime longitudinal volatility \
      --m-metadata-file results/regional-ctf/${rgn}-ctf_state-subject-ordination.qza \
      --p-state-column Quarter_num \
      --p-individual-id-column subject_id \
      --p-default-group-column NERR \
      --p-default-metric PC1 \
      --o-visualization results/regional-ctf/${rgn}-ctf_state_subject_ordination_longitudinal-volatility.qzv
done



find results/regional-ctf -name *_ctf-qurro_plot.qzv -exec cp {} ../jthmiller.github.io/files/results/nerrs/ctf/ \;

find results/regional-ctf -name *_longitudinal-volatility.qzv -exec cp {} ../jthmiller.github.io/files/results/nerrs/feat-volitility \;


    qiime emperor biplot\
      --i-biplot regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/subject_biplot.qza \
      --m-sample-metadata-file  qiime-swmp-corrected-sample-metadata-subject-metadata.tsv \
      --m-feature-metadata-file NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
      --p-number-of-features 20\
      --o-visualization regional_core-diversity/${rgn}_with-repl/${rgn}_gemelli-ctf/${rgn}_subject_biplot.qzv



## Longitudinal volitility
## https://docs.qiime2.org/2020.2/tutorials/longitudinal/

for rgn in SE NE N-Pacific Pacific-Island; do 

  qiime longitudinal feature-volatility \
    --i-table results/${rgn}_genus-table.qza \
    --m-metadata-file metadata.tsv \
    --p-state-column Quarter_num \
    --p-individual-id-column Site_Corrected \
    --p-n-estimators 10 \
    --p-random-state 17 \
    --o-filtered-table results/regional-ctf/${rgn}_filtered-table.qzv \
    --o-feature-importance results/regional-ctf/${rgn}_feature-importance.qzv \
     --o-volatility-plot results/regional-ctf/${rgn}_volatility-plot.qzv \
     --o-accuracy-results results/regional-ctf/${rgn}_accuracy-results.qzv \
     --o-sample-estimator results/regional-ctf/${rgn}_sample-estimator.qzv
done


mv results/regional-ctf/${rgn}_gemelli-ctf/${rgn}_ecam-feat-volatility/volatility_plot.qzv results/feat-volitility/${rgn}_volatility_plot.qzv


cp results/regional-ctf/*_volatility-plot.qzv 

cp results/regional-ctf/*_volatility-plot.qzv ../jthmiller.github.io/files/results/nerrs/feat-volitility/



## Exclude Pac Island
qiime feature-table filter-samples \
  --i-table results/NERRS_18s_euks_hum_samples-table.qza  \
  --m-metadata-file metadata.tsv \
  --p-where "Region!='Pacific-Island'" \
  --o-filtered-table results/NO-island_NERRS_18s-table.qza

qiime gemelli ctf\
    --i-table results/NO-island_NERRS_18s-table.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --p-state-column Quarter_num\
    --p-individual-id-column Site_Corrected\
    --output-dir results/regional-ctf/NO-island_gemelli-ctf

qiime qurro loading-plot \
    --i-ranks results/regional-ctf/NO-island_gemelli-ctf/subject_biplot.qza \
    --i-table results/NO-island_NERRS_18s-table.qza \
    --m-sample-metadata-file metadata.tsv\
    --m-feature-metadata-file results/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza \
    --o-visualization results/regional-ctf/NO-island_gemelli-ctf/NO-island_ctf-qurro_plot.qzv

