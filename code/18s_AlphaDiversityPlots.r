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


## Env types
envtype <- factor(get_variable(psk, "NERR"))
## automatic color palette: one color per different sample type  
palette <- hue_pal()(length(levels(envtype)))
## Map sample type to color
tipColor <- col_factor(palette, levels = levels(envtype))(envtype)  ## Change hclust object to phylo object and plot
sample_data(psk)$NERR_Col <- tipColor

sample_data(psk)[,'salinity'] <- NA
ind <- which(sample_data(psk)[,'Sal_Min.'] > 13)
sample_data(psk)[ind,'salinity'] <- 'SW'

ind <- which(sample_data(psk)[,'Sal_Min.'] < 13)
sample_data(psk)[ind,'salinity'] <- 'FW'

ind <- is.na(sample_data(psk)[,'Sal_Min.'])
sample_data(psk)[ind,'salinity'] <- 'unknown'

sample_data(psk)$Quarter <- as.numeric(factor(sample_data(psk)$Quarter_num))


pruned <- function(ps){
  ps1 <- prune_samples(sample_sums(ps) > 500, ps)
  ps1 <- prune_samples(sample_sums(ps1) < 500000, ps1)
  ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
  return(ps1)
}

psk <- pruned(psk)
psv <- pruned(psv)





    metric = 'Shannon'
    prn <- sample_names(psk)[sample_data(psk)[,'NERR']==NERR]
    psk_sub <- prune_samples(prn, psk)
    psk_sub <- merge_samples(psk_sub, "NERR_SITE_QTR")
    # repair variables that were damaged during merge (coerced to numeric)
    sample_data(psk_sub)$NERR_SITE_QTR <- factor(sample_names(psk_sub))
    sample_data(psk_sub)$Quarter_num <- as.numeric(factor(sample_data(psk_sub)$Quarter_num))
    sample_data(psk_sub)$Quarter <- factor(sample_data(psk_sub)$Quarter)
    sample_data(psk_sub)$NERR <- NERR
    sample_data(psk_sub)$NERR <- factor(sample_data(psk_sub)$NERR, levels = levels(as.factor(sample_data(psk)$NERR)))
    sample_data(psk_sub)$SITE <- factor(unname(unlist(sapply(sample_names(psk_sub), getvars, colinsub = 'NERR_SITE_QTR', colinmd = 'SITE'))))
    p <- plot_richness(psk_sub, x='Quarter', measures = metric) + 
      xlab('Quarter') + 
      scale_y_continuous(limits=c(0,7.5)) + 
      geom_boxplot(lwd=0.9, alpha=0.7) +
      geom_point( aes(colour = factor(SITE)), size = 4) +
      guides(colour = guide_legend(title = "Sites")) +
      theme(axis.text.x = element_text(angle = 90))