# Reads individual consensus peakStats (MACS2 peak stats for each sample at consensus peaks) files and builds composite matrices (FDR, peak signal, Pvalues).
# Set up for MACS2 narrow peaks but easily extended to process narrow and broad peaks in parallel.

library(ggplot2)
library(reshape2)
library(parallel)
library(dplyr)

rm(list = ls())

filelist.narrow <- list.files(path = "./data/narrowPeaks/peakStats",pattern = "*narrowPeaks.peakStats",full.names = T)
sampsStats <- lapply(setNames(as.list(c("narrow")),c("narrow")), function(pset) {
  print(pset)
  fL <- get(paste("filelist",pset,sep = "."))
  stats <- mclapply(fL, function(f) {
    this.samp <- sub("_PBMC.*","",basename(f))
    print(paste("Now processing ",f," (",which(fL==f)," of ",length(fL),")",sep = ""))
    this.stats <- read.table(f,sep = "\t", quote = "", header = F,na.strings = c(".","-1"),col.names = c("chr","start","end","width","name","score","signal","pvalue","qvalue","olap"))
    this.stats$peakco <- with(this.stats,paste(chr,start,end,sep="_"))
    this.stats$pcover <- 100*(this.stats$olap/this.stats$width)
    this.stats[is.na(this.stats$pcover),]$pcover <- NA
    this.stats[this.stats=="."] <- NA
    this.stats <- this.stats[order(this.stats$peakco,-this.stats$signal),]
    this.stats.simple <- this.stats[!duplicated(this.stats$peakco),c("chr","start","end","peakco","signal","pvalue","qvalue","pcover")]
    colnames(this.stats.simple) <- paste(colnames(this.stats.simple),this.samp,sep = "_")
    colnames(this.stats.simple) <- sub("chr.*","chr",sub("start.*","start",sub("end.*","end",sub("peakco.*","peakco",colnames(this.stats.simple)))))
    return(this.stats.simple)
  },mc.cores = 1)
})
compiled.stats <- mclapply(sampsStats, function(S) Reduce(function(x, y) merge(x, y, by = c("chr","start","end","peakco")), S), mc.cores = 1)
compiled.stats.sorted <- mclapply(compiled.stats, function(S) {
  Y <- S %>% arrange(chr,start,end)
  X <- as.matrix(Y %>% dplyr::select(-chr,-start,-end,-peakco))
  X[is.na(X)] <- 0
  Z <- cbind(Y %>% dplyr::select(chr,start,end,peakco),as.data.frame(X))
  return(Z)
},mc.cores = 1)

compiled.signals <- mclapply(compiled.stats.sorted, function(S) cbind(S %>% select(chr,start,end,peakco),
                                                                      S %>% 
                                                                        select(grep("^signal_",colnames(.))) %>% 
                                                                        mutate(Total.counts=rowSums(. > 0),
                                                                               median_signal=apply(.,1,median,na.rm=T),
                                                                               mean_signal=rowMeans(.),
                                                                               max_signal=apply(.,1,max,na.rm=T),
                                                                               min_signal=apply(.,1,min,na.rm=T))),
                             mc.cores = 1)
compiled.pvalues <- mclapply(compiled.stats.sorted, function(S) cbind(S %>% select(chr,start,end,peakco),
                                                                      S %>% 
                                                                        select(grep("^pvalue_",colnames(.))) %>% 
                                                                        mutate(Total.counts=rowSums(. > 0),
                                                                               median_pvalue=apply(.,1,median,na.rm=T),
                                                                               mean_pvalue=rowMeans(.),
                                                                               max_pvalue=apply(.,1,max,na.rm=T),
                                                                               min_pvalue=apply(.,1,min,na.rm=T))),
                             mc.cores = 1)
compiled.qvalues <- mclapply(compiled.stats.sorted, function(S) cbind(S %>% select(chr,start,end,peakco),
                                                                      S %>% 
                                                                        select(grep("^qvalue_",colnames(.))) %>% 
                                                                        mutate(Total.counts=rowSums(. > 0),
                                                                               median_qvalue=apply(.,1,median,na.rm=T),
                                                                               mean_qvalue=rowMeans(.),
                                                                               max_qvalue=apply(.,1,max,na.rm=T),
                                                                               min_qvalue=apply(.,1,min,na.rm=T))),
                             mc.cores = 1)

## Export

write.table(compiled.signals[["narrow"]],file = paste("./data/narrowPeaks/peakStats/narrow_peakstats_signals.txt",sep = ""),quote = F, sep = "\t", row.names = F, col.names = T)
write.table(compiled.pvalues[["narrow"]],file = paste("./data/narrowPeaks/peakStats/narrow_peakstats_pvalues.txt",sep = ""),quote = F, sep = "\t", row.names = F, col.names = T)
write.table(compiled.qvalues[["narrow"]],file = paste("./data/narrowPeaks/peakStats/narrow_peakstats_qvalues.txt",sep = ""),quote = F, sep = "\t", row.names = F, col.names = T)
save.image(paste("./data/narrowPeaks/peakStats/peakStats.extract.RData",sep=""))


