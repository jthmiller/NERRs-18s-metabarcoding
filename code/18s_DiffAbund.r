https://genoweb.toulouse.inra.fr/~formation/15_FROGS/6-October2016/FROGS_phyloseq_10102016.pdf


library("DESeq2")
packageVersion("DESeq2")

# learn hybrid silva
ps <- qza_to_phyloseq(
  tree="results/NERRS_18s_5_20_24_rooted-tree.qza", 
  features='results/filtered-by-features_NERRS_18s_5_20_24.qza', 
  taxonomy='results/NERRS_18s_CRUX_hybrid_taxonomy.qza',
  metadata='metadata/swmp-sample-metadata.tsv'
  )




## lub the sample data
sample_data(ps)$lub_date <- as.Date(mdy(sample_data(ps)$lub_date))
sample_data(ps)$met <- paste0(tolower(sample_data(ps)$SWMP_ID),'met')
sample_data(ps)$wc <- paste0(tolower(sample_data(ps)$SWMP_ID),'wq')

## attach NERRs metadata
data <- read.csv('metadata/797370.csv', sep=',')
data$StationCode <- gsub(" ",'',data$StationCode)
data$lub_time <- as.Date(mdy_hm(data$DateTimeStamp))



## must be prev in 2 samples above 5 reads
## ps0 <- core(ps, detection = 5, prevalence = 2 / 508, include.lowest = T)
ps1 <- prune_samples(sample_sums(ps) > 500, ps)
ps1 <- prune_samples(sample_sums(ps1) < 500000, ps1)

ps1_euks_tree_01 <- tree_glom(ps1, 0.01)
ps1_euks_tree_01 <- prune_taxa(taxa_sums(ps1_euks_tree_01) > 0, ps1_euks_tree_01)





## taxa only
## ps1_euks <- subset_taxa(ps1, Kingdom=="Eukaryota")

## tax agg
# glom_euks <- tip_glom(ps1_euks, tax_adjust = 1)

# Add a prefix to taxa labels
# ps1.f2 <- format_to_besthit(ps1, prefix = "besthit-")

## Tree based clustering
ps1_euks_tree_01 <- tree_glom(ps1_euks, 0.01)
ps1_euks_tree_01 <- prune_taxa(taxa_sums(ps1_euks_tree_01) > 0, ps1_euks_tree_01)

ps1_euks_tree_025 <- tree_glom(ps1_euks, 0.025)
ps1_euks_tree_025 <- prune_taxa(taxa_sums(ps1_euks_tree_025) > 0, ps1_euks_tree_025)


ps1_euks_tree_025 <- tree_glom(ps1_euks, 0.025)
ps1_euks_tree_025 <- prune_taxa(taxa_sums(ps1_euks_tree_025) > 0, ps1_euks_tree_025)


function()

# Calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm = TRUE){
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}





test_geomean <- function(ps, ps_deSeq, alpha){
    geoMeans <- apply(counts(ps_deSeq), 1, gm_mean)
    ps_deSeq <- estimateSizeFactors(ps_deSeq, geoMeans = geoMeans)
    ps_deSeq <- estimateDispersions(ps_deSeq, fitType = "parametric")
    ps_deSeq <- nbinomWaldTest(ps_deSeq)

    res = results(ps_deSeq, cooksCutoff = FALSE)
    
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}

alpha <- 0.01
region_025 <- phyloseq_to_deseq2(ps1_euks_tree_025, ~ Region)
region_025_results <- test_geomean(ps = ps1_euks_tree_025, ps_deSeq = region_025, alpha = alpha)
head(region_025_results)








require(vegan)

 Constrained Analysis of Principal Coordinates (CAP) tries to :
§  Find associaEons between community composiEon and environmental variables (pH, group)
§  QuanEfy differences between groups of samples
 How it works : Regress a distance matrix against some covariates using the standard R syntax for linear
models.
§  Project distance matrix on metadata variables è communiEes are constrained to depend on metadata
§  Look if constrained distance fit to non constrained distance.
 ## convert sample_data to data.frame

dist.uf <- distance(ps1, method = "unifrac")  ## if not already done  
clust.uf <- as.phylo(hclust(dist.uf, method = "complete"))  


 metadata <- as(sample_data(ps1), "data.frame")  
 
cap <- capscale(dist.uf ~ NERR, data = metadata)
anova <- anova(cap, permutations = 999)

ad <- adonis2(dist.uf ~ NERR*Quarter_num*SITE, data = metadata, perm = 9999)



cap <- capscale(dist.uf ~ NERR*Quarter_num, data = metadata)
anova <- anova(cap, permutations = 999)

cap <- capscale(dist.uf ~ NERR*Quarter_num*SITE, data = metadata)
anova <- anova(cap, permutations = 999)

capscale(dist.uf ~ NERR, data = metadata)


capscale(dist.uf ~ , data = metadata)




cap <- capscale(dist.uf ~ NERR + Quarter_num + SITE + NERR*Quarter_num*SITE, data = metadata)
ad <- adonis2(dist.uf ~ NERR + Quarter_num + SITE + NERR*Quarter_num*SITE, data = metadata, perm = 9999)




plot_net(enterotype, maxdist = 0.3, color = "SeqTech", shape="Enterotype")


ig <- make_network(ps1, dist.fun="bray", max.dist=0.3)

png('~/images/quarter_num_netw.png')
plot_network(ig, ps1, color="Quarter_num", shape="NERR", line_weight=0.4, label=NULL)
dev.off()



p <- plot_richness(food, color = "EnvType", x = "EnvType", measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))

p <- plot_richness(food, color = "EnvType", x = "NERR + ", measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))


p <- p + geom_boxplot(aes(fill = EnvType), alpha=0.2)
plot(p)
