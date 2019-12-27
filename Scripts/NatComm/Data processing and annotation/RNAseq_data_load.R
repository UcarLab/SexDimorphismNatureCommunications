# Preps (loads, filters, annotates) RNAseq base data for downstream analyses

library(edgeR)
library(biomaRt)
library(GenomicRanges)
library(dplyr)

rm(list = ls())

andred <- function(x) {Reduce("&",x)}

statsfile = "filtered_global.stats_RNA.RData"
load(statsfile) # metadata variables for samples used in the study, as separate vectors: samps,sampid,  age.group, sex, age, age.group, season, batch, batchdate, nreps, frailty, etc

# RNAseq data -- data paired with ATACseq
rna.base.pbmc <- read.table("./data/RNAseq/rna_rawcounts.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F)

rna.base.pbmc <- (rna.base.pbmc %>%
                    rename_all(funs(sub("\\_PBMC","",.))))[,c("EnsemblID","GeneName","chr","start","strand","GeneType",sampid)]
libsize.rna.pbmc <- colSums(rna.base.pbmc[,sampid])

# Filter data to include only expressed transcripts (as per criterion set herein)
min.cpm = 1 # minimum CPM to be considered as expressed in a sample
min.expressed = 2 # minimum number of samples that have to show expression consider the gene being expressed sample-wise
expressed.pbmc <- data.frame(1*(cpm(as.matrix(rna.base.pbmc[,sampid]),log = F,lib.size = libsize.rna.pbmc)>=min.cpm)) %>%
  cbind(rna.base.pbmc[,-which(colnames(rna.base.pbmc) %in% sampid)],.) %>%
  mutate(HY.counts=rowSums(.[sampid[age.group=="HY"]]),
         HM.counts=rowSums(.[sampid[age.group=="HM"]]),
         HO.counts=rowSums(.[sampid[age.group=="HO"]]),
         Total.counts=rowSums(.[sampid]))
isExpressed <- expressed.pbmc$Total.counts >= min.expressed

# Updates expression data
rna.pbmc <- rna.base.pbmc %>%
  filter(isExpressed,
         !chr %in% c("chrX","chrY","chrM")
         )
libsize.rna.pbmc <- colSums(rna.pbmc[,sampid])
expressed.pbmc <- expressed.pbmc %>% inner_join(rna.pbmc[,c("EnsemblID","GeneName")],by = c("EnsemblID","GeneName"))

write.table(rna.pbmc,"./data/RNAseq/rna_expressed.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(expressed.pbmc[isExpressed,],"./data/RNAseq/rna_calls.txt",sep = "\t",quote = F,col.names = T,row.names = F)

saveobj <- ls()[grep("pbmc",ls())]
save(list = saveobj,file = "./data/RNAseq/pbmc_RNAseq_full.RData")

