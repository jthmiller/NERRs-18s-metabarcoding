## BARPLOTS for replicates per site 
##################################################################
##################################################################

plot_all_bars <- function(site, phyloseq, outappend){

  subsite <- prune_samples(meta(phyloseq)$Region == site, phyloseq)
  subsite_log = transform_sample_counts(subsite, log)
  subsite_rel <- transform_sample_counts(subsite, function(x) x/sum(x))

  praw <- plot_bar(subsite, "Replicate", fill="Phylum", facet_grid=Site_Corrected~Quarter_num) + theme(text = element_text(size = 12))
  plog <- plot_bar(subsite_log, "Replicate", fill="Phylum", facet_grid=Site_Corrected~Quarter_num) + theme(text = element_text(size = 12))
  prel <- plot_bar(subsite_rel, "Replicate", fill="Phylum", facet_grid=Site_Corrected~Quarter_num) + theme(text = element_text(size = 12))

  outpraw <- paste0(site, paste0(outappend, "_raw-barplot.pdf"))
  outplog <- paste0(site, paste0(outappend, "_log-barplot.pdf"))
  outprel <- paste0(site, paste0(outappend, "_relative-barplot.pdf"))

  ggsave(outpraw, plot=praw, width=15, height=40, units="in")
  ggsave(outplog, plog, width=15, height=40, units="in")
  ggsave(outprel, prel, width=15, height=40, units="in")

}
regions <- unique(meta(pc_fam)$Region) 
lapply(regions, plot_all_bars, phyloseq = pc_fam, outappend = '_crux')
lapply(regions, plot_all_bars, phyloseq = ps_fam, outappend = '_silva')






## not used




#### Barplots #####################################################################
physeq_sitemean <- merge_samples2(learn, "samp_rep", fun_otu = mean)
physeq_sitemean <- transform_sample_counts(physeq_sitemean, floor)
# physeq_sitemean <- prune_samples(sample_sums(physeq_sitemean) > 1000, physeq_sitemean)

physeq_transform  <- transform_sample_counts(physeq_sitemean, function(x) x / sum(x) )
physeq_transform <- prune_taxa(taxa_sums(physeq_transform) > 0, physeq_transform)

plot_bar(physeq_transform, "samp_rep", "Abundance", "Phylum") +
facet_grid(~site, scales="free_x", space = "free_x") +
theme(text = element_text(size = 22))
ggsave("NERRS_18s_relative_barplot.png", width=30, height=20, units="in")

plot_bar(physeq_sitemean, "samp_rep", "Abundance", "Phylum") +
facet_grid(~site, scales="free_x", space = "free_x") +
theme(text = element_text(size = 18))
## ggsave("NERRS_18s_absolute_barplot.pdf", width=20, height=20, units="in")


ggsave("NERRS_18s_top5_relative_barplot.png", width=12, height=12, units="in")

### most abundant phyla barplots
phylum.sum = tapply(taxa_sums(physeq_sitemean), tax_table(physeq_sitemean)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
physeq_sitemean_top5 = prune_taxa((tax_table(physeq_sitemean)[, "Phylum"] %in% top5phyla), physeq_sitemean)

plot_bar(physeq_sitemean_top5, "samp_rep", "Abundance", "Phylum") +
facet_grid(~site, scales="free_x", space = "free_x")
ggsave("NERRS_18s_top5_absolute_barplot.pdf", width=20, height=20, units="in")

physeq_sitemean_top5_transform  <- transform_sample_counts(physeq_sitemean_top5, function(x) x / sum(x) )
physeq_sitemean_top5_transform <- prune_taxa(taxa_sums(physeq_sitemean_top5_transform) > 0, physeq_sitemean_top5_transform)

plot_bar(physeq_sitemean_top5_transform, "samp_rep", "Abundance", "Phylum") +
facet_grid(~site, scales="free_x", space = "free_x") +
theme(text = element_text(size = 18))
ggsave("NERRS_18s_top5_relative_barplot.png", width=12, height=12, units="in")


