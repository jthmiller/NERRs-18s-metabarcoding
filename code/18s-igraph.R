

## testing
site <- regions[1]
phyloseq <- pcrux
outappend <- '_crux'
outdir <- 'images/network-plots/'


plot_all_ig <- function(site, phyloseq, outappend, outdir){

  print(site)
  region <- prune_samples(meta(phyloseq)$Region == site, phyloseq)
  #region_log = transform_sample_counts(region, log)
  #region_rel <- transform_sample_counts(region, function(x) x/sum(x))

  region_genus <- tax_glom(region, "Genus")
  region_genus <- merge_samples2(region_genus, "NERR_SITE_QTR", fun_otu = mean)
  region_genus <- transform_sample_counts(region_genus, floor)

  region <- merge_samples2(region, "NERR_SITE_QTR", fun_otu = mean)
  region <- transform_sample_counts(region, floor)

  ig_asv_unifrac <- make_network(region, max.dist=0.65, distance='unifrac')
  ig_asv_jacard <- make_network(region, max.dist=0.95, distance='jaccard')
  ig_genus_unifrac <- make_network(region_genus, max.dist=0.65, distance='unifrac')
  
  #ig_genus <- make_network(region, max.dist=0.45, distance='dpcoa')
  #igraph_options(verbose = TRUE)

  ig_asv_unifrac <- plot_network(ig_asv_unifrac, region, color="salinity", shape="NERR", curved=TRUE)
  ig_asv_jacard <- plot_network(ig_asv_jacard, region, color="salinity", shape="NERR", curved=TRUE)
  ig_genus_unifrac <- plot_network(ig_genus_unifrac, region, color="salinity", shape="NERR", curved=TRUE)

  out_asv_unifrac <- paste0(site, paste0(outappend, "asv_unifrac_ig-network.pdf"))
  out_asv_jacard <- paste0(site, paste0(outappend, "asv_jacard_ig-network.pdf"))
  out_genus_unifrac <- paste0(site, paste0(outappend, "genus_unifrac_ig-network.pdf"))

  out_asv_unifrac <- paste0(outdir,out_asv_unifrac)
  out_asv_jacard <- paste0(outdir,out_asv_jacard)
  out_genus_unifrac <- paste0(outdir,out_genus_unifrac)


  ggsave(out_asv_unifrac, ig_asv_unifrac, width=10, height=10, units="in")
  ggsave(out_asv_jacard, ig_asv_jacard, width=10, height=10, units="in")
  ggsave(out_genus_unifrac, ig_genus_unifrac, width=10, height=10, units="in")
  print(site)

}
regions <- unique(meta(pcrux)$Region) 
lapply(regions, plot_all_ig, phyloseq = pcrux, outappend = '_crux_', outdir = 'images/network-plots/')


lapply(regions, plot_all_ig, phyloseq = psilva, outappend = '_silva_')


