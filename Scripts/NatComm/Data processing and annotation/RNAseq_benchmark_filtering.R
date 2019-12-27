# Finds bechmark peaks for sample filtering, filters data, emits diagnostic plots

rm(list = ls())

benchmark.transcripts <- function(reads=NULL,
                                  expressed=NULL,
                                  anno=NULL,
                                  samps=NULL,
                                  libsize=NULL,
                                  group=NULL,
                                  max.logFC.quantile=0.1,
                                  min.maxct.quantile=0.9,
                                  min.readcount.quantile=0.9,
                                  filter.thresh=0.95) {
  require(edgeR)
  require(GenomicRanges)
  require(ggplot2)
  require(plyr)
  require(dplyr)
  andred <- function(x) {Reduce("&",x)}
  
  print("Invariant (benchmark) transcripts defined as those transcripts:")
  print(paste("---with logFC (HO-HY) in the bottom ",max.logFC.quantile*100,"th percentile of the distribution",sep = ""))
  print(paste("---with normalized counts in the top ",(1-min.maxct.quantile)*100,"th percentile of the distribution",sep = ""))
  print(paste("---with read counts over samples in the top ",(1-min.readcount.quantile)*100,"th percentile of the distribution (i.e., trasncripts called in at least ",floor(min.readcount.quantile*length(group))," out of ",length(group)," samples)",sep = ""))
  print(paste("Filter threshold is set to exclude samples for which ",filter.thresh*100,"% of the benchmark transcripts fail to be called.",sep = ""))
  
  reads.tmm <- setNames(as.data.frame(as.matrix(reads[,samps]) %*% diag(calcNormFactors(as.matrix(reads[,samps]),lib.size = libsize))),samps)
  reads.cpm <-cpm(reads.tmm,lib.size = libsize,log = F)
  benchfilters <- data.frame(merge(expressed[,samps],data.frame(
    mean.HY=apply(reads.cpm[,samps[which(group=="HY")]],1,mean),
    mean.HM=apply(reads.cpm[,samps[which(group=="HM")]],1,mean),
    mean.HO=apply(reads.cpm[,samps[which(group=="HO")]],1,mean),
    maxct=apply(reads.cpm,1,max),
    read.cv=apply(reads.cpm,1,function(x) sd(x)/mean(x)),
    read.count=apply(expressed[,samps],1,sum),
    row.names = rownames(reads.cpm)
  ),by = "row.names"),row.names = 1)
  benchfilters$logFC <- with(benchfilters,log2(1+mean.HO)-log2(1+mean.HY))
  benchfilters <- benchfilters[rownames(reads.cpm),]
  is.invariant <- apply(cbind(
    abs(benchfilters$logFC)<=quantile(abs(benchfilters$logFC),max.logFC.quantile),
    benchfilters$maxct>=quantile(benchfilters$maxct,min.maxct.quantile),
    benchfilters$read.count>=floor(min.readcount.quantile*length(group))
  ),1,andred)
  print(paste("For the chosen parameters, a total of",sum(is.invariant),"invariant (benchmark) reads were defined."))
  reads.invariant <- reads.cpm[is.invariant,]
  benchreads <- data.frame(merge(anno,reads.invariant,by.x = "EnsemblID",by.y = "row.names"),row.names = 1)
  bench.counts <- colSums(expressed[rownames(reads.invariant),samps])
  bench.prop <- bench.counts/sum(is.invariant)
  bench.filter <- bench.prop>filter.thresh
  print(paste("A total of",sum(bench.filter),"out of",length(bench.filter),"samples passed benchmark."))
  if (sum(!bench.filter)>0) print(paste("Failing samples:",paste(samps[!bench.filter],collapse = " ")))
  
  # Export benchmark results
  saveobj.1 <- c("benchreads","bench.counts","bench.prop","reads.invariant","bench.filter")
  save(list = saveobj.1,file = "./data/RNAseq/RNA_full_benchmark.analysis.RData")
  write.table(benchreads,"./data/RNAseq/benchmark.transcripts_RNA_full.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  print(paste("Benchmark transcript IDs written to","./data/RNAseq/benchmark.transcripts_RNA.txt"))
  
  # Export Filtered data
  age <- age[bench.filter]
  age.group <- droplevels(age.group[bench.filter])
  batch <- batch[bench.filter]
  batchdate <- droplevels(batchdate[bench.filter])
  cell.type <- cell.type[bench.filter]
  cmv.status <- cmv.status[bench.filter]
  cmv.titer <- cmv.titer[bench.filter]
  depth <- depth[bench.filter]
  sampid <- names(bench.filter[bench.filter])
  samps <- samps[bench.filter]
  season <- season[bench.filter]
  sex <- droplevels(sex[bench.filter])
  cell.composition <- cell.composition[,bench.filter]
  
  saveobj.2 <- c("age","age.group","batch","batchdate","cell.type","cmv.status","cmv.titer","depth","sampid","samps","season","sex","cell.composition")
  save(list = saveobj.2,file = "filtered_global.stats_RNA_filtered.RData")
  
  # Remove peaks that have no calls after removing samples
  expressed.pbmc <- expressed %>% 
    select(which(!colnames(expressed) %in% names(bench.filter[!bench.filter]))) %>%
    mutate(EnsemblID=rownames(.),
           HY.counts=rowSums(.[sampid[age.group=="HY"]]),
           HM.counts=rowSums(.[sampid[age.group=="HM"]]),
           HO.counts=rowSums(.[sampid[age.group=="HO"]]),
           Total.counts=rowSums(.[sampid])) %>%
    data.frame(.,row.names = "EnsemblID",stringsAsFactors = F)
  rna.pbmc <- reads %>% 
    select(-which(colnames(reads) %in% names(bench.filter[!bench.filter])))
  libsize.rna.pbmc <- colSums(rna.pbmc[,sampid])
  anno.rna.pbmc <- anno %>% inner_join(data.frame(EnsemblID=rownames(expressed.pbmc),stringsAsFactors=F),by="EnsemblID")
  
  saveobj.3 <- c("expressed.pbmc","rna.pbmc","libsize.rna.pbmc","anno.rna.pbmc")
  fdata_file = "./data/RNAseq/pbmc_RNA_full_filtered.RData"
  save(list = saveobj.3,file = fdata_file)
  print(paste("Filtered transcript data exported to",fdata_file))
  
  bench.out <- list(counts=bench.counts,props=bench.prop)
  return(bench.out)
}

### Run the function
load("filtered_global.stats_RNA.RData")
load("./data/RNAseq/pbmc_RNAseq_full.RData") # output from RNAseq_data_load.R

# Restrict metadata to available samples
sampid <- intersect(sampid,sub("_PBMC","",colnames(rna.base.pbmc)))
samps = sampid
age.group <- age.group[sampid]
age <- age[sampid]
sex <- sex[sampid]
season <- season[sampid]
batch <-batch[sampid]
batchdate <- batchdate[sampid]
cell.composition <- cell.composition[,sampid]

anno <- rna.pbmc[,!colnames(rna.pbmc) %in% sampid]
expressed.pbmc <- data.frame(expressed.pbmc[,c("EnsemblID",sampid)],row.names = "EnsemblID")
rna.pbmc <- data.frame(rna.pbmc[,c("EnsemblID",sampid)],row.names = "EnsemblID")
libsize.rna.pbmc <- colSums(rna.pbmc)

rna.bench <- benchmark.transcripts(reads=rna.pbmc,
                                   expressed=expressed.pbmc,
                                   anno=anno,
                                   samps=sampid,
                                   libsize=libsize.rna.pbmc,
                                   group=age.group,
                                   max.logFC.quantile=0.25,
                                   min.maxct.quantile=0.9,
                                   min.readcount.quantile=0.9,
                                   filter.thresh=0.925)
