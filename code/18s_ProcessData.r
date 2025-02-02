## process physeq 

require(qiime2R)
require(phyloseq)
require(tidyverse)
require(speedyseq)
require(plyr)
require(ggpubr)
require(microbiome)
require(microbiomeutilities)
require(RColorBrewer)
require(lubridate)

## what about clustering, then idedntifying the singleton taxa, and remove those asvs

# learn hybrid silva

# learn hybrid silva
psilva <- qza_to_phyloseq(
  tree="qiime-files/input/NERRS_18s_9_12_24_rooted-tree.qza", 
  features='qiime-files/input/NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='qiime-files/input/NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza',
  metadata='qiime-files/input/fixed-qiime-swmp-corrected-sample-metadata.tsv'
  )

pcrux <- qza_to_phyloseq(
  tree="qiime-files/input/NERRS_18s_9_12_24_rooted-tree.qza", 
  features='qiime-files/input/NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='qiime-files/input/NERRS_18s_vsearch_taxonomy_10accepts_90perc_CRUX.qza',
  metadata='qiime-files/input/fixed-qiime-swmp-corrected-sample-metadata.tsv'
  )


tax_table(psilva)[,'Kingdom'] <- gsub('d__','', tax_table(psilva)[,'Kingdom'])
tax_table(pcrux)[,'Kingdom'] <- gsub('d__','', tax_table(pcrux)[,'Kingdom'])

pcrux = subset_taxa(pcrux, Kingdom=="Eukaryota")
psilva = subset_taxa(psilva, Kingdom=="Eukaryota")

pcrux <- subset_taxa(pcrux, (Species!="Homo sapiens") | is.na(Species))
psilva <- subset_taxa(psilva, (Species!="Homo_sapiens") | is.na(Species))

## Drop very low samples to prevent missing data from driving the signal.
## Dropped any sample that has less than 50 asvs
psilva <- prune_samples(sample_sums(psilva)>=50, psilva)
pcrux <- prune_samples(sample_sums(pcrux)>=50, pcrux)

# Aggregate at the family level 
pc_fam <- tax_glom(pcrux, "Family")
ps_fam <- tax_glom(psilva, "Family")

#pc0 <- tip_glom(pcrux, 0.1, tax_adjust = 0)
pc1 <- tip_glom(pcrux, 0.05, tax_adjust = 1)
#pc2 <- tip_glom(pcrux, 0.1, tax_adjust = 2)

pcphy = tax_glom(pcrux, "Phylum")
pctax <- tax_glom(pcrux, "Genus")
pctree <- tree_glom(pcrux, resolution = 0.05, criterion = "max_tip_depth", tax_adjust = 1L)
pctip <- tip_glom(pcrux, 0.05, tax_adjust = 1)

#filter_taxa(pctree5, function(x) sum(x > 1) > 1, TRUE)
#filter_taxa(pctree5, function(x) sum(x > 1) > 1, TRUE)
#pctree_none <- tree_glom(pcrux, resolution = 0.05, criterion = "max_tip_depth", tax_adjust = 0L)
## rowSums(otu_table(pcrux) > 0) > 1

## filter taxa only found in a single sample
pcrux <- filter_taxa(pcrux, function(x) sum(x > 1) > 1, TRUE)
psilva <- filter_taxa(psilva, function(x) sum(x > 1) > 1, TRUE)

## After dropping taxa, reroot the tree with ape 
phy_tree(psilva) <- ape::multi2di(phy_tree(psilva))
phy_tree(pcrux) <- ape::multi2di(phy_tree(pcrux))




otucnt <- rowSums(otu_table(pcrux) > 0)
pdf('otu-hist.pdf')
hist(otucnt[otucnt>0], breaks = 500)
dev.off()