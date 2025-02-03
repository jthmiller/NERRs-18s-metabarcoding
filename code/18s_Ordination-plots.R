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
require(ape)
require(scales)

## cd /home/users/jtm1171/NERRs/18s
psk <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_CRUX_hybrid_taxonomy.qza',
  metadata='phylo-swmp-sample-metadata.tsv'
  )

psv <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_vsearch_taxonomy_10accepts_90perc.qza',
  metadata='phylo-swmp-sample-metadata.tsv'
  )

### Metadata
sample_data(psk)$Quarter <- as.numeric(factor(sample_data(psk)$Quarter_num))


## AW updated salinity 
##sample_data(psk)[ which(  sample_data(psk)[,'Sal_Min.'] > 40 ),'Sal_Min.'] <- NA
##sample_data(psk)[,'salinity'] <- NA
##ind <- which(sample_data(psk)[,'Sal_Min.'] > 13)
##sample_data(psk)[ind,'salinity'] <- 'SW'
##ind <- which(sample_data(psk)[,'Sal_Min.'] < 13)
##sample_data(psk)[ind,'salinity'] <- 'FW'


## Env types
envtype <- factor(get_variable(psk, "NERR"))
## automatic color palette: one color per different sample type  
palette <- hue_pal()(length(levels(envtype)))
## Map sample type to color
tipColor <- col_factor(palette, levels = levels(envtype))(envtype)  ## Change hclust object to phylo object and plot
sample_data(psk)$NERR_Col <- tipColor



pruned <- function(ps){
  ps1 <- prune_samples(sample_sums(ps) > 500, ps)
  ps1 <- prune_samples(sample_sums(ps1) < 500000, ps1)
  ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
  return(ps1)
}

psk <- pruned(psk)
psv <- pruned(psv)

############################################################################
############################################################################
## taxa only
psk_euks <- subset_taxa(psk, Kingdom=="Eukaryota")

### Regional
psk_euks_tf <- transform_sample_counts(psk_euks, function(x) x/sum(x))
psk_euks_tf <- prune_taxa(taxa_sums(psk_euks_tf) > 0, psk_euks_tf)

psk_euks_tf_Atlantic <- subset_samples(psk_euks_tf, Ocean=='Atlantic')     
psk_euks_tf_Gulf <- subset_samples(psk_euks_tf, Ocean=='Gulf')  
psk_euks_tf_Pacific <- subset_samples(psk_euks_tf, Ocean=="Pacific") 
### Regional



#ps1_euks_tree_01 <- tree_glom(ps1, 0.01)
#ps1_euks_tree_01 <- prune_taxa(taxa_sums(ps1_euks_tree_01) > 0, ps1_euks_tree_01)







## tax agg
# glom_euks <- tip_glom(ps1_euks, tax_adjust = 1)

# Add a prefix to taxa labels
# ps1.f2 <- format_to_besthit(ps1, prefix = "besthit-")

## Tree based clustering
psk_euks_tree_01 <- tree_glom(psk_euks, 0.01)
psk_euks_tree_01 <- prune_taxa(taxa_sums(psk_euks_tree_01) > 0, psk_euks_tree_01)

psk_euks_tf_tree = transform_sample_counts(psk_euks_tree_01, function(x) x/sum(x))
psk_euks_tf_tree <- prune_taxa(taxa_sums(psk_euks_tf_tree) > 0, psk_euks_tf)


#psk_euks_tree_025 <- tree_glom(psk_euks, 0.025)
#psk_euks_tree_025 <- prune_taxa(taxa_sums(psk_euks_tree_025) > 0, psk_euks_tree_025)
#
#psk_euks_tree_025 <- tree_glom(psk_euks, 0.025)
#psk_euks_tree_025 <- prune_taxa(taxa_sums(psk_euks_tree_025) > 0, psk_euks_tree_025)


### Ordination #########################################

## Transform to even sampling depth.
## physeq_tf = transform_sample_counts(physeq, function(x) 1E6 * x/sum(x))
psk_euks_tf = transform_sample_counts(psk_euks, function(x) x/sum(x))
psk_euks_tf <- prune_taxa(taxa_sums(psk_euks_tf) > 0, psk_euks_tf)

### Requires taxa to be present in at least 2 samples
psk_euks_tf3 <- filter_taxa(psk_euks, function(x) sum(x > 1) > 1, TRUE)
psk_euks_tf3 = transform_sample_counts(psk_euks_tf3, function(x) x/sum(x))
psk_euks_tf3 <- prune_taxa(taxa_sums(psk_euks_tf3) > 0, psk_euks_tf3)
# physeq <- transform_sample_counts(physeq, floor)

#### ALL SITES
ordu = ordinate(psk_euks_tf, "PCoA", "unifrac", weighted=TRUE)
wns <- plot_ordination(psk_euks_tf, ordu, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')
ggsave("~/NERRS_18s_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

ordu = ordinate(psk_euks_tf, "PCoA", "unifrac", weighted=TRUE)
wns <- plot_ordination(psk_euks_tf, ordu, color="North_South", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')
ggsave("~/NERRS_18s_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

ordu = ordinate(psk_euks_tf, "PCoA", "unifrac", weighted=F)
uwns <- plot_ordination(psk_euks_tf, ordu, color="North_South", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')
ggsave("~/NERRS_18s_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")

ordu = ordinate(psk_euks_tf, "PCoA", "unifrac", weighted=TRUE)
w <- plot_ordination(psk_euks_tf, ordu, color="Ocean", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

ordu = ordinate(psk_euks_tf, "PCoA", "unifrac", weighted=F)
uw <- plot_ordination(psk_euks_tf, ordu, color="Ocean", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("~/NERRS_18s_PCoA_unifrac_500reads_min.pdf", width=10, height=10, units="in")
#### ALL SITES









### Regional PCoA_unifrac_weighted
ordu_Atlantic = ordinate(psk_euks_Atlantic, "PCoA", "unifrac", weighted=TRUE)
wns_Atlantic <- plot_ordination(psk_euks_Atlantic, ordu_Atlantic, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Weighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

psk_euks_Gulf_tf <- transform_sample_counts(psk_euks_Gulf, function(x) x/sum(x))
psk_euks_Gulf_tf <- prune_taxa(taxa_sums(psk_euks_Gulf_tf) > 0, psk_euks_Gulf_tf)
ordu_Gulf = ordinate(psk_euks_Gulf_tf, "PCoA", "unifrac", weighted=TRUE)
wns_Gulf <- plot_ordination(psk_euks_Gulf_tf, ordu_Gulf, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Weighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Gulf_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

psk_euks_Pacific_tf <- transform_sample_counts(psk_euks_Pacific, function(x) x/sum(x))
psk_euks_Pacific_tf <- prune_taxa(taxa_sums(psk_euks_Pacific_tf) > 0, psk_euks_Pacific_tf)
ordu_Pacific = ordinate(psk_euks_Pacific_tf, "PCoA", "unifrac", weighted=TRUE)
wns_Pacific <- plot_ordination(psk_euks_Pacific_tf, ordu_Pacific, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Weighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Pacific_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")
### Regional PCoA_unifrac_weighted


### Regional PCoA_unifrac_unweighted
ordu_Atlantic_uw <- ordinate(psk_euks_Atlantic_tf, "PCoA", "unifrac", weighted=FALSE)
wns_Atlantic <- plot_ordination(psk_euks_Atlantic_tf, ordu_Atlantic_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Unweighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")

ordu_Atlantic_uw <- ordinate(psk_euks_Atlantic_tf, "PCoA", "unifrac", weighted=FALSE)
wns_Atlantic <- plot_ordination(psk_euks_Atlantic_tf, ordu_Atlantic_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Unweighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")

ordu_Atlantic_uw <- ordinate(psk_euks_Atlantic_tf, "PCoA", "unifrac", weighted=FALSE)
wns_Atlantic <- plot_ordination(psk_euks_Atlantic_tf, ordu_Atlantic_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Unweighted PCoA, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")
### Regional PCoA_unifrac_unweighted







### tree vs ASV
## ASV
ordu_Atlantic = ordinate(psk_euks_tf_Atlantic, "PCoA", "unifrac", weighted=TRUE)
wns_Atlantic <- plot_ordination(psk_euks_tf_Atlantic, ordu_Atlantic, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Weighted PCoA, Unifrac, ASVs')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

## Tree based clustering
psk_euks_tf_Atlantic_tree_01 <- tree_glom(psk_euks_tf_Atlantic, 0.01)
ordu_Atlantic_tree = ordinate(psk_euks_tf_Atlantic_tree_01, "PCoA", "unifrac", weighted=TRUE)
wns_Atlantic_tree <- plot_ordination(psk_euks_tf_Atlantic_tree_01, ordu_Atlantic_tree, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Weighted PCoA Tree, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")


## ASV
ordu_Pacific = ordinate(psk_euks_tf_Pacific, "PCoA", "unifrac", weighted=TRUE)
wns_Pacific <- plot_ordination(psk_euks_tf_Pacific, ordu_Pacific, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Weighted PCoA, Unifrac, ASVs')
ggsave("~/NERRS_18s_Pacific_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

## Tree based clustering
psk_euks_tf_Pacific_tree_01 <- tree_glom(psk_euks_tf_Pacific, 0.01)
ordu_Pacific_tree = ordinate(psk_euks_tf_Pacific_tree_01, "PCoA", "unifrac", weighted=TRUE)
wns_Pacific_tree <- plot_ordination(psk_euks_tf_Pacific_tree_01, ordu_Pacific_tree, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Weighted PCoA Tree, Unifrac')
ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

## ASV
ordu_Gulf = ordinate(psk_euks_tf_Gulf, "PCoA", "unifrac", weighted=TRUE)
wns_Gulf <- plot_ordination(psk_euks_tf_Gulf, ordu_Gulf, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Weighted PCoA, Unifrac, ASVs')
#ggsave("~/NERRS_18s_Gulf_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

## Tree based clustering
psk_euks_tf_Gulf_tree_01 <- tree_glom(psk_euks_tf_Gulf, 0.01)
ordu_Gulf_tree = ordinate(psk_euks_tf_Gulf_tree_01, "PCoA", "unifrac", weighted=TRUE)
wns_Gulf_tree <- plot_ordination(psk_euks_tf_Gulf_tree_01, ordu_Gulf_tree, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Weighted PCoA Tree, Unifrac')
#ggsave("~/NERRS_18s_Atlantic_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

ggarrange(wns_Atlantic, wns_Atlantic_tree, wns_Pacific, wns_Pacific_tree, wns_Gulf, wns_Gulf_tree, ncol = 3, nrow = 3)
ggsave("~/NERRS_18s_PCoA_unifrac_500reads_min.pdf", width=10, height=10, units="in")
### tree vs ASV




### ASV weighted v unweighted PCOA regional
ordu_Atlantic = ordinate(psk_euks_tf_Atlantic, "PCoA", "unifrac", weighted=TRUE)
ordu_Atlantic_uw = ordinate(psk_euks_tf_Atlantic, "PCoA", "unifrac", weighted=FALSE)
ordu_Pacific = ordinate(psk_euks_tf_Pacific, "PCoA", "unifrac", weighted=TRUE)
ordu_Pacific_uw = ordinate(psk_euks_tf_Pacific, "PCoA", "unifrac", weighted=FALSE)
ordu_Gulf = ordinate(psk_euks_tf_Gulf, "PCoA", "unifrac", weighted=TRUE)
ordu_Gulf_uw = ordinate(psk_euks_tf_Gulf, "PCoA", "unifrac", weighted=FALSE)
### ASV weighted v unweighted PCOA regional


## Plot by salinity
wns_Atlantic <- plot_ordination(psk_euks_tf_Atlantic, ordu_Atlantic, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Weighted PCoA, Unifrac, ASVs')

wns_Atlantic_uw <- plot_ordination(psk_euks_tf_Atlantic, ordu_Atlantic_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Unweighted PCoA, Unifrac, ASVs')

wns_Pacific <- plot_ordination(psk_euks_tf_Pacific, ordu_Pacific, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Weighted PCoA, Unifrac, ASVs')

wns_Pacific_uw <- plot_ordination(psk_euks_tf_Pacific, ordu_Pacific_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Unweighted PCoA, Unifrac, ASVs')

wns_Gulf <- plot_ordination(psk_euks_tf_Gulf, ordu_Gulf, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Weighted PCoA, Unifrac, ASVs')

wns_Gulf_uw <- plot_ordination(psk_euks_tf_Gulf, ordu_Gulf_uw, color="salinity", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Unweighted PCoA, Unifrac, ASVs')

ggarrange(wns_Atlantic, wns_Atlantic_uw, wns_Pacific, wns_Pacific_uw, wns_Gulf, wns_Gulf_uw, ncol = 2, nrow = 3)
ggsave("~/NERRS_18s_PCoA_unifrac_regional_salinity.pdf", width=10, height=10, units="in")
## Plot by salinity




### Min Salinity Colored
wns_Atlantic <- plot_ordination(psk_euks_tf_Atlantic, ordu_Atlantic, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Weighted PCoA, Unifrac, ASVs')

wns_Atlantic_uw <- plot_ordination(psk_euks_tf_Atlantic, ordu_Atlantic_uw, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Atlantic, Unweighted PCoA, Unifrac, ASVs')

wns_Pacific <- plot_ordination(psk_euks_tf_Pacific, ordu_Pacific, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Weighted PCoA, Unifrac, ASVs')

wns_Pacific_uw <- plot_ordination(psk_euks_tf_Pacific, ordu_Pacific_uw, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Pacific, Unweighted PCoA, Unifrac, ASVs')

wns_Gulf <- plot_ordination(psk_euks_tf_Gulf, ordu_Gulf, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Weighted PCoA, Unifrac, ASVs')

wns_Gulf_uw <- plot_ordination(psk_euks_tf_Gulf, ordu_Gulf_uw, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('Gulf, Unweighted PCoA, Unifrac, ASVs')

ggarrange(wns_Atlantic, wns_Atlantic_uw, wns_Pacific, wns_Pacific_uw, wns_Gulf, wns_Gulf_uw, ncol = 2, nrow = 3)
ggsave("~/NERRS_18s_PCoA_unifrac_regional_min_sal.pdf", width=10, height=10, units="in")
### Min Salinity Colored