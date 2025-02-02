




## colors
n <- 42
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols <- sample(col_vector, n)



# check healthy
GB_ps1 <- subset_samples(ps1_euks, NERR=="GB")
p_gb <- taxa_distribution(GB_ps1, color.taxa = cols) + 
  theme_biome_utils() + 
  labs(title = "Great Bay")

# check CRC
WE_ps1 <- subset_samples(ps1_euks, NERR=="WE")
p_we <- taxa_distribution(WE_ps1, color.taxa = cols) + 
  theme_biome_utils() + 
  labs(title = "Wells")

png('~/phylum-dist.png',width = 1360)
p_gb / p_we + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A")
dev.off()




## 'Cyanobacteria/Chloroplast related sequences which can be removed if not expected in the samples'
p <- plot_read_distribution(ps, groups = "NERR", plot.type = "density")

myTaxa = names(sort(taxa_sums(ps1_euks), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, GlobalPatterns)
plot(phy_tree(ex1), show.node.label = TRUE)

myTaxa = names(sort(rowSums(otu_table(ps1_euks)>0), decreasing = TRUE)[1:10])
ex1 = prune_taxa(myTaxa, ps1_euks)
png('~/tree_top_phyla.png')
plot(phy_tree(ex1), show.node.label = TRUE)
dev.off()

png('~/tree_top_phyla.png',width = 1360, height = 1000)
plot_tree(ex1, color = "NERR", label.tips = "Class", ladderize = "left", justify = "left" , size = "Abundance")
dev.off()

png('~/tree_top_phyla.png',width = 1360, height = 1000)
plot_tree(ex1, color = "NERR", nodelabf=nodeplotblank, label.tips = "Class", ladderize = "left", size = "Abundance")
dev.off()
plot.margin = 0.2
"OTU"



## Tree of verts only filtered
verts = subset_taxa(ps1_euks, Phylum=="Vertebrata")
verts = format_to_besthit(verts, prefix = NULL)
tax_table(verts)[,'best_hit'] <- gsub(".*:","",rownames(otu_table(verts)))
speedyseq::plot_tree(verts, label.tips = "Species",size="abundance", color="site", text.size =4,nodelabf=nodeplotboot(80,0,3)) 
ggsave("~/NERRS_18s_filt_verts.pdf", width=10, height=10, units="in")
## + coord_polar(theta="y")


## dominant taxa
p0.gen <- aggregate_taxa(ps1_euks,"Phylum")
x.d <- dominant_taxa(ps1_euks,level = "Phylum", group="NERR")
head(x.d$dominant_overview, 10)

p0.gen <- aggregate_taxa(ps1_euks,"Family")
x.d <- dominant_taxa(ps1_euks,level = "Family", group="NERR")
head(x.d$dominant_overview, 10)
