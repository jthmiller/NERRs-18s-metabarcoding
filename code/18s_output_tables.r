
cd ~/NERRs.. 
conda activate qiime2-amplicon-2024.5
export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.422.b05-2.el9.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH'
R 

both composiEonal and phylogeneEc
Unifrac
Weigthed-Unifrac
FracEon of the tree specific to either 1 or 2
FracEon of the diversity specific to 1orto2

- Jaccard higher than Unifrac
    - communiEes taxa are disEnct but phylogeneEcally related


### https://genoweb.toulouse.inra.fr/~formation/15_FROGS/6-October2016/FROGS_phyloseq_10102016.pdf


#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla
## export LD_LIBRARY_PATH='/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH'
## Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH')

Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.422.b05-2.el9.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH')


library(tidyverse)
library(qiime2R)
library(taxize)
require(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(rJava)
library(xlsx)
require(speedyseq)
require(openxlsx)
## Sys.setenv(LD_LIBRARY_PATH = '/usr/lib/jvm/java-1.11.0-openjdk-amd64/lib/server:$LD_LIBRARY_PATH')


### use: this_script.R <feature_table.qza> <taxonomy.qza> <rep-seqs.qza> <ASV_outname.csv> <Species_outname.csv>

options <- commandArgs(trailingOnly = TRUE)


# options <- c("results/NERRS_18s_5_20_24_filtered-table.qza",
#   "results/NERRS_18s_5_20_24_hybrid_taxonomy.qza",
#   "results/NERRS_18s_5_20_24_rep-seqs.qza",
#   "results/NERRS_18s_5_20_24",
#   "metadata/18s_sample_metadata_5-20.tsv")


 #options <- c(
 #  "results/filtered_HIDAR-COI_table_wControls.qza",
 #  "results/HIDAR-COI_midori_hybrid_taxonomy.qza",
 #  "results/HIDAR-COI_rep-seqs.qza",
 #  "results/HDAR_COI_midori"
 #  )

##options <- c(
##  'qiime_out/HD24-1-18SNX041524_table.qza',
##   'qiime_out/HD24-1-18SNX041524_06042024_hybrid_taxonomy.qza',
##    'qiime_out/HD24-1-18SNX041524_rep-seqs.qza',
##    'qiime_out/HD24-1-18SNX041524_06042024',
##    'qiime_out/HD24-1-18SNX041524_dns_export/metadata.tsv'
##)

#options <- c(
#  'results/MBON_TERNS_JUN17_table_filtered.qza',
#  'results/MBON_TERNS_JUN17_hybrid_taxonomy.qza',
#  'results/MBON_TERNS_JUN17_rep-seqs.qza',
#  'MBON_TERNS_JUN17',
#  'metadata/MBON_TERNS_JUN17_fitler-metadata.tsv')
#
#options <- c(
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_table.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_06212024_hybrid_taxonomy.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_rep-seqs.qza',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624',
#  'runs/HD24-1-CO1NX041624/qiime_out/HD24-1-CO1NX041624_dns_export/metadata.tsv') 
#options <- c(
#  'HD23-RERUN-MFNX051424_table.qza',
#  'HD23-RERUN-MFNX051424_07162024_hybrid_taxonomy.qza',
#  'HD23-RERUN-MFNX051424_rep-seqs.qza',
#  'HD23-RERUN-MFNX051424',
#  'HD23-RERUN-MFNX051424_dns_export/metadata.tsv')
#

options <- c(
  'NERRS_18s_9_12_24_filtered-table.qza',
  'NERRS_18s_CRUX_hybrid_taxonomy.qza',
  'NERRS_18s_5_20_24_rep-seqs.qza',
  'NERRS_18s_CRUX',
  'phylo-swmp-sample-metadata.tsv',
  'NERRS_18s_5_20_24_rooted-tree.qza'
  )

## plot sequences per sample
physeq <- qza_to_phyloseq(
  features = options[1], 
  taxonomy = options[2],
  tree = options[6]
  )

physeq <- qza_to_phyloseq(
  features = options[1], 
  taxonomy = options[2],
  metadata = 'phylo-swmp-sample-metadata.tsv'
  )


reps <- read_qza(options[3])$data
tax <- read_qza(options[2])$data
feats <- read_qza(options[1])$data
rownames(tax) <- tax$Feature.ID

## d__ doesn't work
tax_table(physeq)[,'Kingdom'] <- gsub("d__",'',tax_table(physeq)[,'Kingdom'])
tax_table(physeq)[,'Kingdom'] <- gsub("tax=",'',tax_table(physeq)[,'Kingdom'])


### ASV table
ASVs_out <- data.frame(
  ASV_ID = rownames(tax_table(physeq)),
  tax_table(physeq),
  Confidence = tax[rownames(tax_table(physeq)),'Confidence'],
  ASV_Sequence = as.character(reps[rownames(tax_table(physeq))]),
  num_samples_pos = rowSums(otu_table(physeq)>0),
  otu_table(physeq),
  check.names = FALSE
  )

############################################################################################
# Species table
## Collapse taxa (should not have sequences- multiple sequences per taxa)
otus <- tax_glom(physeq, taxrank=rank_names(physeq)[7], NArm=F, bad_empty=c(NA, "", " ", "\t"))

otus_out <- data.frame(
  ASV_ID = rownames(tax_table(otus)),
  tax_table(otus),
  num_samples_pos = rowSums(otu_table(otus)>0),
  total_reads = rowSums(otu_table(otus)),
  otu_table(otus),
  check.names = FALSE
  )
otus_out <- otus_out[,-1]

# read count
#otu_reads_out <- data.frame(cbind(colnames(otu_table(physeq)), colSums(otu_table(physeq))))
############################################################################################


############################################################################################
# Table on tree based clustering
tree_otus <- tree_glom(physeq, 0.025)
tree_otus <- prune_taxa(taxa_sums(tree_otus) > 0, tree_otus)

tree_out <- data.frame(
  ASV_ID = rownames(tax_table(tree_otus)),
  tax_table(tree_otus),
  num_samples_pos = rowSums(otu_table(tree_otus)>0),
  total_reads = rowSums(otu_table(tree_otus)),
  otu_table(tree_otus),
  check.names = FALSE
  )
tree_out <- tree_out[,-1]




############################################################################################
# Table on clutsering tips
tip_otus_agnes <- tip_glom(physeq,  h = 0.025, hcfun = cluster::agnes)
tip_otus_agnes <- prune_taxa(taxa_sums(tip_otus_agnes) > 0, tip_otus_agnes)

tip_otus_hclust <- tip_glom(physeq,  h = 0.025, hcfun = hclust)
tip_otus_hclust <- prune_taxa(taxa_sums(tip_otus_hclust) > 0, tip_otus_hclust)


tip_out_agnes <- data.frame(
  ASV_ID = rownames(tax_table(tip_otus_agnes)),
  tax_table(tip_otus_agnes),
  num_samples_pos = rowSums(otu_table(tip_otus_agnes)>0),
  total_reads = rowSums(otu_table(tip_otus_agnes)),
  otu_table(tip_otus_agnes),
  check.names = FALSE
  )
tip_out_agnes <- tip_out_agnes[,-1]


tip_out_hclust <- data.frame(
  ASV_ID = rownames(tax_table(tip_otus_hclust)),
  tax_table(tip_otus_hclust),
  num_samples_pos = rowSums(otu_table(tip_otus_hclust)>0),
  total_reads = rowSums(otu_table(tip_otus_hclust)),
  otu_table(tip_otus_hclust),
  check.names = FALSE
  )
tip_out_hclust <- tip_out_hclust[,-1]

############################################################################################


# attach metadata
metadata_out <- read.table(options[5], sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = "@")


ASV_file <- paste0(options[4],"_ASV_table.csv")
otu_file <- paste0(options[4],"_classifier-based-otus_table.csv")
tree_file <- paste0(options[4],"_tree-based-otus_table.csv")
tip_file <- paste0(options[4],"_tip-based-otus_table.csv")

ASV_file <- paste0('~/',ASV_file)
otu_file <- paste0('~/',otu_file)
tree_file <- paste0('~/',tree_file)
tip_file <- paste0('~/',tip_file)


write.csv(ASVs_out, file = ASV_file, row.names=FALSE, quote=FALSE)
write.csv(otus_out, file = otu_file, row.names=FALSE, quote=FALSE)
write.csv(tree_out, file = tree_file, row.names=FALSE, quote=FALSE)
write.csv(tip_out, file = tip_file, row.names=FALSE, quote=FALSE)

xlsx_out <- paste0(options[4],'.xlsx')
xlsx_out <- paste0('~/',xlsx_out)

out <- list(ASVs_out, otus_out, tree_out, tip_out_agnes, metadata_out)
names(out) <- c("ASVs","taxonomy-based-otus","tree-based-otus","clustering-based-otus","metadata")
openxlsx::write.xlsx(out, file = xlsx_out, rowNames=FALSE)




# write.xlsx(ASVs_out, file = xlsx_out, sheetName="ASVs", row.names=FALSE)
# write.xlsx(otus_out, file = xlsx_out, sheetName="OTUs", append=TRUE, row.names=FALSE)
# write.xlsx(reads_out, file = xlsx_out, sheetName="readcount", append=TRUE, row.names=FALSE)
# write.xlsx(tax_out, file = xlsx_out, sheetName="taxonomy", append=TRUE, row.names=FALSE)

# openxlsx::write.xlsx(otus_out, file = xlsx_out, sheetName="OTUs", append=TRUE, rowNames=FALSE)
# openxlsx::write.xlsx(reads_out, file = xlsx_out, sheetName="readcount", append=TRUE, rowNames=FALSE)
# openxlsx::write.xlsx(tax_out, file = xlsx_out, sheetName="taxonomy", append=TRUE, rowNames=FALSE)