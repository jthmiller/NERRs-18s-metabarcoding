





ps <- gb_all_tf
function(ps,NERR,feat)

samps <- rownames(sample_data(ps))[which(sample_data(ps)$NERR == NERR)]

otu_abun <- otu_table(ps)[feat,samps]





sample <- rownames(sample_data(ps))[1]


## ## NERRs 18s Final Analysis
## conda config --prepend pkgs_dirs /home/users/jtm1171/.conda/pkgs
## conda config --prepend envs_dirs /home/users/jtm1171/.conda/envs
## 
## ##### NEW ENVIRONMENT FOR Qiime2R
## conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
## conda activate qiime2-amplicon-2024.5
## conda install bioconda::bioconductor-phyloseq
## conda install r-qiime2r
## conda install bioconda::bioconductor-phyloseq
## conda install conda-forge::r-tidyverse
## conda install conda-forge::r-devtools


# ## Install in R
# R
# remotes::install_github("adeverse/ade4")
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
# install.packages("remotes")
# remotes::install_github("mikemc/speedyseq")
# install.packages("ggpubr")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("microbiome")
# devtools::install_github("microsud/microbiomeutilities")
# install.packages('patchwork')

# results/NERRS_18s_5_20_24_rep-seqs.qza
# results/NERRS_18s_5_20_24_table.qza
# results/NERRS_18s_5_20_24_hybrid_taxonomy.qza
# results/NERRS_18s_5_20_24_aligned-rep-seqs.qza
# results/NERRS_18s_5_20_24_masked-aligned-rep-seqs.qza
# results/NERRS_18s_5_20_24_unrooted-tree.qza
# results/NERRS_18s_5_20_24_rooted-tree.qza
# results/NERRS_18s_5_20_24_filtered-table.qza
# results/18s_sample_metadata_5-20.tsv_phyloseq

# cd /home/users/jtm1171/old-home/watts/data/NERR/18s
# conda activate qiime2-amplicon-2024.5
# R 

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


# learn hybrid crux
ps <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_vsearch_taxonomy_10accepts_90perc_CRUX.qza',
  metadata='../../old-home/watts/data/NERR/18s/metadata/NERRs_18s_Sample_metadata.tsv'
  )

## lub the sample data
sample_data(ps)$lub_date <- as.Date(mdy(sample_data(ps)$lub_date))
sample_data(ps)$met <- paste0(tolower(sample_data(ps)$SWMP_ID),'met')
sample_data(ps)$wc <- paste0(tolower(sample_data(ps)$SWMP_ID),'wq')
sample_data(ps)$nut <- paste0(tolower(sample_data(ps)$SWMP_ID),'nut')

### swmp data wc and met
data <- read.csv('../../old-home/watts/data/NERR/18s/metadata/797370.csv', sep=',')
data$StationCode <- gsub(" ",'',data$StationCode)
data$lub_time <- as.Date(mdy_hm(data$DateTimeStamp))
wq_data <- data[ which(data$StationCode %in% sample_data(ps)$wc),]
met_data <- data[ which(data$StationCode %in% sample_data(ps)$met),]

data <- read.csv('../../old-home/watts/data/NERR/18s/metadata/500389.csv', sep=',')
data$StationCode <- gsub(" ",'',data$StationCode)
data$lub_time <- as.Date(mdy_hm(data$DateTimeStamp))
nut_data <- data[ which(data$StationCode %in% sample_data(ps)$nut),]

rm(data)


wq <- c('Temp','SpCond','Sal','DO_pct','DO_mgl','Depth','pH','Turb')
nut <- c('PO4F','NH4F','NO23F','CHLA_N')
met <- c('ATemp','WSpd','MaxWSpd','TotPAR','TotPrcp')

## remove bad data


summary(nut_data[,nut])


summary(wq_data[,wq])
wq_data[,'Temp'][which(wq_data[,'Temp'] < -25)] <- NA
wq_data[,'SpCond'][which(wq_data[,'SpCond'] < 0)] <- NA
wq_data[,'SpCond'][which(wq_data[,'SpCond'] > 100)] <- NA
wq_data[,'Sal'][which(wq_data[,'Sal'] > 50)] <- NA
wq_data[,'DO_pct'][which(wq_data[,'DO_pct'] > 120)] <- NA
wq_data[,'DO_pct'][which(wq_data[,'DO_pct'] < 0)] <- NA
wq_data[,'DO_mgl'][which(wq_data[,'DO_mgl'] > 50)] <- NA
wq_data[,'DO_mgl'][which(wq_data[,'DO_mgl'] < -5)] <- NA
wq_data[,'Depth'][which(wq_data[,'Depth'] < 0)] <- NA
wq_data[,'pH'][which(wq_data[,'pH'] > 20)] <- NA
wq_data[,'pH'][which(wq_data[,'pH'] < 0)] <- NA
wq_data[,'Turb'][which(wq_data[,'Turb'] < 0)] <- NA

summary(met_data[,met])
met_data[,'TotPAR'][which(met_data[,'TotPAR'] < 0)] <- NA
# hist(met_data[,'TotPrcp'][which(met_data[,'TotPrcp'] > 0)])


#swmps <- unique(tolower(sample_data(ps)$SWMP_ID))
#swmps <- swmps[!swmps == '']
#swmpdata_codes <- 
#matchswmp <- lapply(swmps, grep, x=unique(data$StationCode),value=T)
#names(matchswmp) <- swmps
#nodat <- which(unlist(lapply(lapply(samp_day_data,dim),'[[',1)) == 0)



get_day_data <- function(sample, ps, data, col){
  dat <- sample_data(ps)
  ind <- which(rownames(dat) == sample)
  swmp_site <- as.character(dat[ind,col])
  sampdate <- as.Date(dat$lub_date[ind])
  var <- paste0('days_since_',col)

  if (length( which(data$StationCode == swmp_site)) == 0){ 
    out <- data.frame(t(rep(NA, times= length(colnames(data)))))
    colnames(out) <- colnames(data)
    diff <- NA
  } else {
    out <- data[which(data$StationCode == swmp_site),]
    time_vec <- abs(sampdate - out$lub_time)
    
    if (any(as.numeric(time_vec) < 14)){
      diff <- as.numeric(min(time_vec))
      out <- out[ which(time_vec == diff),]
    } else {
      out <- data.frame(t(rep(NA, times= length(colnames(data)))))
      colnames(out) <- colnames(data)
      diff <- NA
    }
  }
  outnames <- c(colnames(out),var)
  out <- data.frame(out,diff)
  colnames(out) <- outnames
  return(out)

}
 
samps <- rownames(sample_data(ps))
wc_samp_day_data <- lapply(samps, get_day_data, ps = ps, data = wq_data, col = 'wc')
nut_samp_day_data <- lapply(samps, get_day_data, ps = ps, data = nut_data, col = 'nut')
met_samp_day_data <- lapply(samps, get_day_data, ps = ps, data = met_data, col = 'met')

names(wc_samp_day_data) <- samps
names(nut_samp_day_data) <- samps
names(met_samp_day_data) <- samps



### Get avg to attach for simple plot and testing
wc <- c('Temp','SpCond','Sal','DO_pct','DO_mgl','Depth','pH','Turb')
nut <- c('PO4F','NH4F','NO23F','CHLA_N')
met <- c('ATemp','WSpd','MaxWSpd','TotPAR','TotPrcp')

stat <- c("Min.","Median","Mean","Max.")



sumize <- function(samp ,data, params, stats){
  out <- lapply(params, function(param, sam = samp){ 
    setNames(summary(sam[,param])[stats],paste(param,stats, sep='_')) 
    })
  days_since <- grep('days_since_',colnames(samp), value = T)
  return(setNames(c(unlist(out),samp[1,days_since]),c(names(unlist(out)),days_since))) 
}
nut_summary <- lapply(nut_samp_day_data, sumize, params = nut, stats = stat)
met_summary <- lapply(met_samp_day_data, sumize, params = met, stats = stat)
wc_summary <- lapply(wc_samp_day_data, sumize, params = wc, stats = stat)

sum_param <- data.frame(do.call(rbind,wc_summary),do.call(rbind,met_summary),do.call(rbind,nut_summary))
rownames(sum_param) <- names(met_samp_day_data)




### Write new metadata
sample_data(ps)[,'salinity'] <- NA
ind <- which(sample_data(ps)[,'Sal_Min.'] > 13)
sample_data(ps)[ind,'salinity'] <- 'SW'

ind <- which(sample_data(ps)[,'Sal_Min.'] < 13)
sample_data(ps)[ind,'salinity'] <- 'FW'

ind <- is.na(sample_data(ps)[,'Sal_Min.'])
sample_data(ps)[ind,'salinity'] <- 'unknown'


### write the metadata
df <- data.frame(cbind(sample_data(ps), sum_param[rownames(sample_data(ps)),]))
SampleID <- rownames(df)
df <- cbind(SampleID,df)
## fix qiime error
colnames(df)[which(colnames(df) == "sample_name")] <- 'prev_sample_name'
write.table(df,'swmp-corrected-sample-metadata.tsv',sep = '\t', quote=F, row.names=F)
### write the metadata



# learn hybrid silva
ps <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza',
  metadata='phylo-qiime-swmp-corrected-sample-metadata.tsv'
  )


Quarter_num




### write the metadata
df <- data.frame(cbind(sample_data(ps), sum_param[rownames(sample_data(ps)),]))

df <- data.frame(sample_data(ps))
SampleID <- rownames(df)
df <- cbind(SampleID,df)
df[,'Quarter_num'] <- as.numeric(factor(as.numeric(df[,'Quarter_num'])))
## fix qiime error
colnames(df)[which(colnames(df) == "sample_name")] <- 'prev_sample_name'
write.table(df,'swmp-corrected-sample-metadata.tsv',sep = '\t', quote=F, row.names=F)
### write the metadata


write.table(tax_table(ps), 'silva-18s-tax.csv',sep = ',', quote=F, row.names=T)


repseqs <- read_qza('NERRS_18s_5_20_24_rep-seqs.qza')$data
as.character(repseqs['c1816b32f8c0f0caa7caf17cbe7a98e0'])


as.character(repseqs['6e887399a2514e2224f8414fc0e789d4'])






# learn hybrid silva
ps <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza',
  metadata='phylo-qiime-swmp-corrected-sample-metadata.tsv'
  )

repseqs <- read_qza('NERRS_18s_5_20_24_rep-seqs.qza')$data

out <- as.character(repseqs[rownames(otu_table(ps))])
out <- data.frame(cbind(rownames(out),out))
write.table(out,'sequences.tsv',sep = ',', quote=F, row.names=T)

head(read_qza('NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza')$data)








# learn hybrid silva
ps <- qza_to_phyloseq(
  tree="NERRS_18s_9_12_24_rooted-tree.qza", 
  features='NERRS_18s_9_12_24_filtered-table.qza', 
  taxonomy='NERRS_18s_vsearch_taxonomy_10accepts_90perc-silva.qza',
  metadata='phylo-qiime-swmp-corrected-sample-metadata.tsv'
  )

### SOME SAMPLE NAME PROBLEMS STILL 
#SSSE, all have extra reps 
#SFCC, Q1 has extra reps 
#all SS sites missing for quarter 1. 

### Drop these (sampling usually takes place at low tide. Thgese are high):
## sample_data(psilva)[ which(sample_data(psilva)[,'Site_Corrected']=='SFCC' & sample_data(psilva)[,'Quarter_num']==1),]
sample_data(ps)[ 'SFCCHw0530231','Site_Corrected'] <- 'SFCCH'
sample_data(ps)[ 'SFCCHw0530232','Site_Corrected'] <- 'SFCCH'
sample_data(ps)[ 'SFCCHw0530233','Site_Corrected'] <- 'SFCCH'

#sample_data(ps)[ which(sample_data(ps)[,'Site_Corrected']=='SFCC'),]

sample_data(pcrux)[ which(sample_data(pcrux)[,'Site_Corrected']=='SSSE'),]

### SSSA, SSWS both changed to SSSE (SSBH SSCH SSSE  SSVA  SSWI)
## sample_data(ps)[ which(sample_data(ps)[,'Site_Corrected']=='SSSE'),]

#SSBH Boat house 9
#SSCH Charls Bridge 9
#SSSE Sengstacken Arm *18*
#SSVA Valino Island 9
#SSWI Winchester Arm 9
#
#SSSE is made up from:
#SSSA (Q2-Q4) is probably SSSE? 
#SSWS (Q2-Q4) maybe SSWI only different?

sample_data(ps)[ 'SSWSw0222241','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw0222242','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw0222243','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw1128231','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw1128232','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw1128233','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw0823231','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw0823232','Site_Corrected'] <- 'SSWS'
sample_data(ps)[ 'SSWSw0823233','Site_Corrected'] <- 'SSWS'

sample_data(ps)[which( meta(ps)$Region == 'SE'),'Region'] <- 'Gulf'

md <- meta(ps)
sample_data(ps)[ ,'Quarter_num'] <- as.numeric(factor(as.numeric(md$Quarter_num)))
sample_data(ps)[ ,'NERR_SITE_QTR'] <- paste(md$Site_Corrected, paste0('Q',  md$Quarter_num), sep="_")


### write the metadata
df <- meta(ps)
SampleID <- rownames(df)
df <- cbind(SampleID,df)
write.table(df,'fixed-qiime-swmp-corrected-sample-metadata.tsv',sep = '\t', quote=F, row.names=F)