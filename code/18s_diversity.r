
cd /home/users/jtm1171/old-home/watts/data/NERR/18s
conda activate qiime2-amplicon-2024.5
R 

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
# learn hybrid silva

ps <- qza_to_phyloseq(
  tree="results/NERRS_18s_5_20_24_rooted-tree.qza", 
  features='results/filtered-by-features_NERRS_18s_5_20_24.qza', 
  taxonomy='results/NERRS_18s_CRUX_hybrid_taxonomy.qza',
  metadata='metadata/swmp-sample-metadata.tsv'
  )


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

pruned <- function(ps){
  ps1 <- prune_samples(sample_sums(ps) > 500, ps)
  ps1 <- prune_samples(sample_sums(ps1) < 500000, ps1)
  ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
  return(ps1)
}

psk <- pruned(psk)
psv <- pruned(psv)



psk_tree <- tree_glom(psk, 0.025)
psk_tree <- prune_taxa(taxa_sums(psk_tree) > 0, psk_tree)

psk


ord.bray <- ordinate(psk, method = "MDS", distance = "bray")
ord.jacard <- ordinate(psk, method = "MDS", distance = "jaccard")
ord.unifrac <- ordinate(psk, method = "MDS", distance = "unifrac")
ord.wunifrac <- ordinate(psk, method = "MDS", distance = "wunifrac")

 plot(hclust(dist.uf, method = "single"))

source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R" )

p <- plot_ordination(ps1, ord.bray, color = "Region")
 ## add title and plain background  
 
p <- p + theme_bw() + ggtitle("MDS + BC")  




dist.bc <- distance(ps1, method = "unifrac") 
png('~/images/unifrac_clust_single.png', width = 4000)
 plot(hclust(dist.bc, method = "single"))
dev.off()




## Env types
envtype <- factor(get_variable(ps1, "NERR"))
## automatic color palette: one color per different sample type  
palette <- hue_pal()(length(levels(envtype)))
## Map sample type to color
tipColor <- col_factor(palette, levels = levels(envtype))(envtype)  ## Change hclust object to phylo object and plot



dist.uf <- distance(ps1, method = "unifrac")  ## if not already done  
clust.uf <- as.phylo(hclust(dist.uf, method = "complete"))  

png('~/images/complete_color_clust.png', width = 5000)
par(mar = c(0, 0, 2, 0))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "complete")
dev.off()

clust.uf.single <- as.phylo(hclust(dist.uf, method = "single"))  
png('~/images/single_color_clust.png', width = 5000)
par(mar = c(0, 0, 2, 0))
plot(clust.uf.single, tip.color = tipColor, direction = "downwards", main = "single")
dev.off()

clust.uf.ward <- as.phylo(hclust(dist.uf, method = "ward.D2"))  
png('~/images/ward.D2_color_clust.png', width = 5000)
par(mar = c(0, 0, 2, 0))
plot(clust.uf.ward, tip.color = tipColor, direction = "downwards", main = "single")
dev.off()


dist.bc <- distance(ps, method = "bray")




richness <- estimate_richness(ps1_euks_tree_01)


png('~/all_richness.png')
plot_richness(ps1_euks_tree_01)
dev.off()


physeq_sub <- prune_samples(sample_data(ps1_euks_tree_01)$Quarter_num == 5, ps1_euks_tree_01)

png('~/Q1_shannon_richness.png')
plot_richness(physeq_sub, x="NERR", measures="Shannon", color = "NERR")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
dev.off()



physeq_sub <- prune_samples(sample_data(ps1_euks_tree_01)$Quarter_num == 8, ps1_euks_tree_01)

png('~/Q4_shannon_richness.png')
plot_richness(physeq_sub, x="NERR", measures="Shannon", color = "NERR")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
dev.off()



######


############# DIVERSITY ##############
############# DIVERSITY ##############
## ps1_euks: High/Low samples filtered only




### Ordination #########################################

## Transform to even sampling depth.
## physeq_tf = transform_sample_counts(physeq, function(x) 1E6 * x/sum(x))
ps1_euks_tf = transform_sample_counts(ps1_euks, function(x) x/sum(x))
ps1_euks_tf <- prune_taxa(taxa_sums(ps1_euks_tf) > 0, ps1_euks_tf)

### Requires taxa to be present in at least 2 samples
ps1_euks_tf3 <- filter_taxa(ps1_euks, function(x) sum(x > 1) > 1, TRUE)
ps1_euks_tf3 = transform_sample_counts(ps1_euks_tf3, function(x) x/sum(x))
ps1_euks_tf3 <- prune_taxa(taxa_sums(ps1_euks_tf3) > 0, ps1_euks_tf3)



ps1_euks_tf = transform_sample_counts(ps1_euks_tree_01, function(x) x/sum(x))
ps1_euks_tf <- prune_taxa(taxa_sums(ps1_euks_tf) > 0, ps1_euks_tf)

# physeq <- transform_sample_counts(physeq, floor)

## PCoA unifrac
#ordu = ordinate(physeq_tf, "PCoA", "unifrac", weighted=TRUE)
#plot_ordination(physeq_tf, ordu, color="samp_rep", shape="site") +
#scale_shape_manual(values=seq(0,9))
#ggsave("NERRS_18s_PCoA_unifrac_300reads_3sampTaxa.pdf", width=20, height=20, units="in")

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=TRUE)
wns <- plot_ordination(ps1_euks_tf, ordu, color="North_South", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')
ggsave("~/NERRS_18s_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=F)
uwns <- plot_ordination(ps1_euks_tf, ordu, color="North_South", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')
ggsave("~/NERRS_18s_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")


ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=TRUE)
w <- plot_ordination(ps1_euks_tf, ordu, color="Ocean", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=F)
uw <- plot_ordination(ps1_euks_tf, ordu, color="Ocean", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("~/NERRS_18s_PCoA_unifrac_500reads_min.pdf", width=10, height=10, units="in")








ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=TRUE)
wns <- plot_ordination(ps1_euks_tf, ordu, color="Quarter_num", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')
# ggsave("~/NERRS_18s_PCoA_unifrac_weighted.pdf", width=20, height=20, units="in")

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=F)
uwns <- plot_ordination(ps1_euks_tf, ordu, color="Quarter_num", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')
#ggsave("~/NERRS_18s_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=TRUE)
w <- plot_ordination(ps1_euks_tf, ordu, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

ordu = ordinate(ps1_euks_tf, "PCoA", "unifrac", weighted=F)
uw <- plot_ordination(ps1_euks_tf, ordu, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("NERRS_18s_PCoA_unifrac_500reads_min_Quarter_Region.pdf", width=10, height=10, units="in")

### unweighted, not normalized

ps1_euks_tree_01

ordu = ordinate(ps1_euks_tree_01, "PCoA", "unifrac", weighted=F)
uwns <- plot_ordination(ps1_euks_tree_01, ordu, color="North_South", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ordu = ordinate(ps1_euks_tree_01, "PCoA", "unifrac", weighted=F)
uwoc <- plot_ordination(ps1_euks_tree_01, ordu, color="Ocean", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ordu = ordinate(ps1_euks_tree_01, "PCoA", "unifrac", weighted=F)
uwqn <- plot_ordination(ps1_euks_tree_01, ordu, color="Quarter_num", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')
#ggsave("~/NERRS_18s_PCoA_unifrac_unweighted.pdf", width=20, height=20, units="in")

ordu = ordinate(ps1_euks_tree_01, "PCoA", "unifrac", weighted=F)
uwre <- plot_ordination(ps1_euks_tree_01, ordu, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(uwns,uwoc,uwqn,uwre, ncol = 2, nrow = 2)
ggsave("~/NERRS_18s_PCoA_unifrac_500reads_unweighted.pdf", width=10, height=10, units="in")










##### BARS ###################################
 p <- plot_bar(ps1, fill = "Phylum")
 ## add facets
 p <- p + facet_wrap(~NERR, scales = "free_x", nrow = 1)
 
 png('~/images/barplot.png',width=3000)
 plot(p)
 dev.off()





p <- plot_composition(ps1,sample.sort="Quarter_num", otu.sort="abundance", x.label= 'SITE', numberOfTaxa = 30, fill = "Phylum", group_by='NERR')

p <- plot_composition(ps1,sample.sort="Quarter_num", otu.sort="abundance", x.label= 'SITE', numberOfTaxa = 30, fill = "Phylum")


 ## plot facetting
 p <- p + facet_wrap(~NERR, scales = "free_x", nrow = 10)  
 
 p <- p + facet_grid(.~NERR+SITE)

 png('~/images/bar-plot-ints.png')
 plot(p)
dev.off()
 #############################################





data.table::melt(sample_data(ps1))

mdf <- melt(df,id.vars=rownames(df))


ggplot(mdf, aes( x=Quarter_num, y=value, colour=variable, group=variable )) + 
  geom_line() +
  scale_color_manual(values=c("y1"="black","y2"="red","y3"="orange")) +
  scale_linetype_manual(values=c("y1"="solid","y2"="solid","y3"="dashed"))



means <- grep('_Mean', colnames(sample_data(ps1)), value=T)
df <- data.frame(sample_data(ps1))
df <- df[ ,c("sample_name","Quarter_num","NERR","NERR_SITE_QTR", means)]



df_lf <- 


df %>% pivot_longer(!sample_name, names_to = "variable", values_to = 'value')



    cols = sample_names(physeq_besthit), 
    names_to = 'sample-id', 
    values_to = 'Found')





png('~/temp_mean.png', width = 750, height = 200)
ggplot( aes(  x = as.factor(Quarter_num), y = Temp_Mean), data = df) + geom_boxplot() + facet_grid(.~NERR) + scale_y_continuous(limits = c(0, 35))
dev.off()


png('~/pH_Mean.png', width = 750, height = 200)
ggplot( aes(  x = as.factor(Quarter_num), y = pH_Mean), data = df) + geom_boxplot() + facet_grid(.~NERR) + scale_y_continuous(limits = c(6,9))
dev.off()


png('~/Sal_Mean.png', width = 750, height = 200)
ggplot( aes(  x = as.factor(Quarter_num), y = Sal_Mean), data = df) + geom_boxplot() + facet_grid(.~NERR) + scale_y_continuous(limits = c(10, 40))
dev.off()

png('~/NO23F_Mean.png', width = 750, height = 200)
ggplot( aes(  x = as.factor(Quarter_num), y = NO23F_Mean), data = df) + geom_boxplot() + facet_grid(.~NERR) + scale_y_continuous(limits = c(0.002, 1.4))
dev.off()








metadata <- read_q2metadata("sample-metadata.tsv")
shannon <- read_qza("shannon_vector.qza")

shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged


metadata<-
metadata %>% 
  left_join(shannon)

metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`days-since-experiment-start`, y=shannon, color=`body-site`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Shannon Diversity") +
  theme_q2r() + # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="Body Site") # use different color scale which is color blind friendly
  ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches








ps1




## no normalization, no clustering
orduW = ordinate(ps1, "PCoA", "unifrac", weighted=TRUE)
orduUW = ordinate(ps1, "PCoA", "unifrac", weighted=F)


wns <- plot_ordination(ps1, orduW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwns <- plot_ordination(ps1, orduUW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

w <- plot_ordination(ps1, orduW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uw <- plot_ordination(ps1, orduUW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("NERRS_18s_PCoA_unifrac_500reads_min_TEMP_Region_ASVs_RAW.pdf", width=10, height=10, units="in")






## no normalization, clustering
orduW = ordinate(ps1_tree, "PCoA", "unifrac", weighted=TRUE)
orduUW = ordinate(ps1_tree, "PCoA", "unifrac", weighted=F)

wns <- plot_ordination(ps1_tree, orduW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwns <- plot_ordination(ps1_tree, orduUW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

w <- plot_ordination(ps1_tree, orduW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uw <- plot_ordination(ps1_tree, orduUW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("~/images/NERRS_18s_PCoA_unifrac_500reads_min_TEMP_Region_ASVs_CLUSTERED.pdf", width=10, height=10, units="in")









## no normalization, no clustering, cold quarters only
ps1_cold <- subset_samples(ps1, Quarter_num > 6)
sample_data(ps1_cold)[ which(sample_data(ps1_cold)[,'Temp_Median'] < -20) ,'Temp_Median'] <- NA
sample_data(ps1_cold)[which(sample_data(ps1_cold)[,'Temp_Median'] < -20),'Temp_Median']  <- NA


orduW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=TRUE)
orduUW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=F)

wns <- plot_ordination(ps1_cold, orduW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwns <- plot_ordination(ps1_cold, orduUW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

w <- plot_ordination(ps1_cold, orduW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uw <- plot_ordination(ps1_cold, orduUW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("~/images/NERRS_18s_PCoA_unifrac_500reads_min_TEMP_Region_ASVs_COLD.pdf", width=10, height=10, units="in")





## no normalization, clustering, cold quarters only
ps1_cold <- subset_samples(ps1_tree, Quarter_num > 6)
sample_data(ps1_cold)[ which(sample_data(ps1_cold)[,'Temp_Median'] < -20) ,'Temp_Median'] <- NA

orduW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=TRUE)
orduUW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=F)

wns <- plot_ordination(ps1_cold, orduW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwns <- plot_ordination(ps1_cold, orduUW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

w <- plot_ordination(ps1_cold, orduW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uw <- plot_ordination(ps1_cold, orduUW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

ggarrange(wns,uwns,w,uw, ncol = 2, nrow = 2)
ggsave("~/images/NERRS_18s_PCoA_unifrac_500reads_min_TEMP_Region_CLUSTERS_COLD.pdf", width=10, height=10, units="in")





### shannon diversity over quarters
p <- plot_richness(ps1_tree, x='Quarter_num', measures="Shannon", color = "NERR")+
  geom_boxplot(aes(group = Quarter_num), alpha=0.6)+ 
  geom_point(aes(group = Quarter_num), alpha=0.6)+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) + 
  scale_x_discrete() +
  facet_wrap(~NERR, nrow = 2)

#### Shannon Diversity
png('~/images/quarters_shannon_richness.png', width = 1000, height = 400)
p
dev.off()









## no normalization, no clustering, cold quarters only
ps1_cold <- subset_samples(ps1, Quarter_num > 6)
sample_data(ps1_cold)[ which(sample_data(ps1_cold)[,'Temp_Median'] < -20) ,'Temp_Median'] <- NA

orduW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=TRUE)
orduUW = ordinate(ps1_cold, "PCoA", "unifrac", weighted=F)

wns <- plot_ordination(ps1_cold, orduW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwns <- plot_ordination(ps1_cold, orduUW, color="Region", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')

w <- plot_ordination(ps1_cold, orduW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uw <- plot_ordination(ps1_cold, orduUW, color="Temp_Median", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')


wl <- plot_ordination(ps1_cold, orduW, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('weighted PCoA unifrac')

uwl <- plot_ordination(ps1_cold, orduUW, color="Sal_Min.", shape="NERR") +
scale_shape_manual(values=seq(0,9)) +
ggtitle('unweighted PCoA unifrac')


ggarrange(wns,uwns,w,uw,wl,uwl, ncol = 2, nrow = 3)
ggsave("~/images/NERRS_18s_PCoA_unifrac_500reads_min_TEMP_Region_ASVs_COLD_7+8.pdf", width=15, height=15, units="in")
