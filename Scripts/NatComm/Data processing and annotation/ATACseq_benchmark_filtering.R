# Finds bechmark peaks for sample filtering, filters data, emits diagnostic plots

rm(list = ls())

benchmark.peaks <- function(peakSet="narrow",
                            peaks=NULL,
                            peakcalls=NULL,
                            anno=NULL,
                            samps=NULL,
                            libsize=NULL,
                            depth=NULL,
                            group=NULL,
                            max.logFC.quantile=0.1,
                            min.maxct.quantile=0.9,
                            min.peakcount.quantile=0.9,
                            filter.thresh=0.95,
                            show.plots=TRUE) {
  require(edgeR)
  require(GenomicRanges)
  require(ggplot2)
  require(plyr)
  require(dplyr)
  andred <- function(x) {Reduce("&",x)}
  
  peakSet = peakSet[1]
  print(paste("Peak set being benchmarked:",toupper(peakSet)))
  print("Invariant (benchmark) peaks defined as those peaks:")
  print(paste("---with logFC (HO-HY) in the bottom ",max.logFC.quantile*100,"th percentile of the distribution",sep = ""))
  print(paste("---with normalized counts in the top ",(1-min.maxct.quantile)*100,"th percentile of the distribution",sep = ""))
  print(paste("---with peak counts over samples in the top ",(1-min.peakcount.quantile)*100,"th percentile of the distribution (i.e., peaks called in at least ",floor(min.peakcount.quantile*length(group))," out of ",length(group)," samples)",sep = ""))
  print(paste("Filter threshold is set to exclude samples for which ",filter.thresh*100,"% of the benchmark peaks fail to be called.",sep = ""))
  
  peaks.tmm <- setNames(as.data.frame(as.matrix(peaks[,samps]) %*% diag(calcNormFactors(as.matrix(peaks[,samps]),lib.size = libsize))),samps)
  rownames(peaks.tmm) <- peaks$peakco
  benchfilters <- data.frame(merge(peakcalls[,samps],data.frame(
    mean.HY=apply(peaks.tmm[,samps[which(group=="HY")]],1,mean),
    mean.HO=apply(peaks.tmm[,samps[which(group=="HO")]],1,mean),
    maxct=apply(peaks.tmm,1,max),
    peak.cv=apply(peaks.tmm,1,function(x) sd(x)/mean(x)),
    peak.count=apply(peakcalls[,samps],1,sum),
    row.names = rownames(peaks.tmm)
  ),by = "row.names"),row.names = 1)
  benchfilters$logFC <- with(benchfilters,log2(mean.HO)-log2(mean.HY))
  benchfilters <- benchfilters[rownames(peaks.tmm),]
  is.invariant <- apply(cbind(
    abs(benchfilters$logFC)<quantile(abs(benchfilters$logFC),max.logFC.quantile),
    benchfilters$maxct>quantile(benchfilters$maxct,min.maxct.quantile),
    benchfilters$peak.count>=floor(min.peakcount.quantile*length(group))
  ),1,andred)
  print(paste("For the chosen parameters, a total of",sum(is.invariant),"invariant (benchmark) peaks were defined."))
  peaks.invariant <- peaks.tmm[is.invariant,]
  benchpeaks <- as.data.frame(sort(GRanges(data.frame(merge(anno,peaks.invariant,by.x = "peakco",by.y = "row.names"),row.names = 1))))
  benchpeaks <- setNames(benchpeaks,sub("seqnames","chr",colnames(benchpeaks)))[,!colnames(benchpeaks)=="strand"]
  bench.counts <- colSums(peakcalls[rownames(peaks.invariant),samps])
  bench.prop <- bench.counts/sum(is.invariant)
  bench.filter <- bench.prop>filter.thresh
  print(paste("A total of",sum(bench.filter),"out of",length(bench.filter),"samples passed benchmark."))
  if (sum(!bench.filter)>0) print(paste("Failing samples:",paste(samps[!bench.filter],collapse = " ")))
  
  # Export benchmark results
  saveobj.1 <- c("benchpeaks","bench.counts","bench.prop","peaks.invariant","bench.filter")
  save(list = saveobj.1,file = paste("./data/",peakSet,"Peaks/",peakSet,"_benchmark.analysis.RData",sep = ""))
  write.table(benchpeaks,paste("./data/",peakSet,"Peaks/benchmark.peaks_",peakSet,".txt",sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
  print(paste("Benchmark peak locations written to",paste("./data/",peakSet,"Peaks/benchmark.peaks_",peakSet,".txt",sep = "")))
  
  if (show.plots) {
    frip <- libsize/depth
    pdf(paste("./data/",peakSet,"Peaks/",peakSet,"_benchmark.filters.pdf",sep = ""))
    bench.plot <- data.frame(frip=frip[samps],bench.prop=bench.prop[samps],peak.set=peakSet,sampid=samps)
    p <- ggplot(aes(frip,bench.prop,label=samps),data=bench.plot) + 
      geom_point(color="white",shape=1,size=2) + 
      geom_text(size=2.5) + 
      ggtitle(paste(peakSet,"peaks")) + 
      ylab("Prop of samples calling benchmark peaks") + 
      xlab("FRiP score") + 
      geom_abline(slope=0,intercept=filter.thresh,color="firebrick1",linetype=2) + 
      geom_abline(slope=0,intercept=filter.thresh-0.1,color="firebrick4",linetype=2)
    print(p)
    
    bench.plot.anno <- rbind.fill(
      data.frame(t(data.frame(prop.table(table(toupper(anno$Annotation))),row.names = 1)),check.names = F),
      data.frame(t(data.frame(prop.table(table(toupper(data.frame(merge(anno,peaks.invariant,by = "row.names"),row.names = 1)$Annotation))),row.names = 1)),check.names = F)
    )
    barplot(as.matrix(bench.plot.anno),beside = T,main = paste(peakSet," peaks (",nrow(peaks.invariant)," benchmark peaks)", sep = ""),col=c("dimgray","tomato"),border=c("dimgray","tomato"),las=2,cex.names=0.6,cex.main = 0.8,cex = 0.7)
    legend(x="topleft",legend=c("All peaks","Benchmark peaks"),fill = c("dimgray","tomato"),border = c("dimgray","tomato"),cex=0.7,bty = "n")
    
    bench.plot.hmm <- rbind.fill(
      data.frame(t(data.frame(prop.table(table(toupper(anno$chromHMMsimple_PBMC))),row.names = 1)),check.names = F),
      data.frame(t(data.frame(prop.table(table(toupper(data.frame(merge(anno,peaks.invariant,by = "row.names"),row.names = 1)$chromHMMsimple_PBMC))),row.names = 1)),check.names = F)
    )
    barplot(as.matrix(bench.plot.hmm),beside = T,main = paste(peakSet," peaks (",nrow(peaks.invariant)," benchmark peaks)", sep = ""),col=c("dimgray","goldenrod1"),border=c("dimgray","goldenrod1"),las=2,cex.names=0.6,cex.main = 0.8,cex = 0.7)
    legend(x="topleft",legend=c("All peaks","Benchmark peaks"),fill = c("dimgray","goldenrod1"),border = c("dimgray","goldenrod1"),cex=0.7,bty = "n")
    dev.off()
    print(paste("Plots saved to ",peakSet,"_benchmark.filters.pdf",sep = ""))
  }
  
  # Export Filtered data
  age <- age[bench.filter]
  age.group <- droplevels(age.group[bench.filter])
  batch <- batch[bench.filter]
  batchdate <- droplevels(batchdate[bench.filter])
  cell.type <- cell.type[bench.filter]
  frailty.index <- frailty.index[bench.filter]
  frail <- frail[bench.filter]
  cmv.status <- cmv.status[bench.filter]
  cmv.titer <- cmv.titer[bench.filter]
  depth <- depth[bench.filter]
  sampid <- names(bench.filter[bench.filter])
  samps <- samps[bench.filter]
  season <- season[bench.filter]
  sex <- droplevels(sex[bench.filter])
  nreps <- nreps[bench.filter]
  cell.composition <- cell.composition[,bench.filter]

  saveobj.2 <- c("age","age.group","batch","cell.type","frailty.index","frail","cmv.status","cmv.titer","depth","sampid","samps","season","sex","batchdate","nreps","cell.composition")
  save(list = saveobj.2,file = paste("filtered_global.stats_",peakSet,"Peaks_filtered.RData",sep = ""))
  
  # Remove peaks that have no calls after removing samples
  peakcalls.pbmc <- peakcalls %>% 
    select(-which(colnames(peakcalls) %in% names(bench.filter[!bench.filter]))) %>%
    mutate(HY.counts=rowSums(.[sampid[age.group=="HY"]]),
           HM.counts=rowSums(.[sampid[age.group=="HM"]]),
           HO.counts=rowSums(.[sampid[age.group=="HO"]]),
           Total.counts=rowSums(.[sampid])) %>%
    filter(Total.counts>=2) %>%
    data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F)
  atac.allpeaks.pbmc <- (peaks %>% 
                           select(-which(colnames(peaks) %in% names(bench.filter[!bench.filter]))) %>%
                           data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F))[rownames(peakcalls.pbmc),]
  libsize.atac.pbmc <- colSums(atac.allpeaks.pbmc[,sampid])
  frips.atac.pbmc <- libsize.atac.pbmc/depth
  atac.allpeaks.anno.pbmc <- anno[rownames(atac.allpeaks.pbmc),]
  atac.allpeaks.hmm.pbmc <- atac.allpeaks.hmm.pbmc[rownames(atac.allpeaks.pbmc),]
  atac.allpeaks.eqtl.pbmc <- atac.allpeaks.eqtl.pbmc %>% filter(peakco %in% atac.allpeaks.anno.pbmc$peakco)
  atac.allpeaks.cspec_eqtl.pbmc <- atac.allpeaks.cspec_eqtl.pbmc %>% filter(peakco %in% atac.allpeaks.anno.pbmc$peakco)
  peakqvalues.pbmc <- peakqvalues.pbmc[rownames(atac.allpeaks.pbmc),] %>%
    select(-which(colnames(peakqvalues.pbmc) %in% names(bench.filter[!bench.filter]))) %>%
    mutate(maxq=apply(.[,sampid],1,max))
  peaksignals.pbmc <- peaksignals.pbmc[rownames(atac.allpeaks.pbmc),] %>%
    select(-which(colnames(peaksignals.pbmc) %in% names(bench.filter[!bench.filter])))
  
  
  saveobj.3 <- c("peakcalls.pbmc","atac.allpeaks.pbmc","libsize.atac.pbmc","frips.atac.pbmc","atac.allpeaks.anno.pbmc","atac.allpeaks.hmm.pbmc","atac.allpeaks.eqtl.pbmc","atac.allpeaks.cspec_eqtl.pbmc","peakqvalues.pbmc","peaksignals.pbmc")
  fdata_file = paste("./data/",peakSet,"Peaks/",peakSet,"Peaks_filtered.RData",sep = "")
  save(list = saveobj.3,file = fdata_file)
  print(paste("Filtered peak data exported to",fdata_file))
  
  bench.out <- list(counts=bench.counts,props=bench.prop)
  return(bench.out)
}

### Run the function
load("filtered_global.stats.RData")
load("./data/narrowPeaks/pbmc_narrowPeaks.RData")
narrow.bench <- benchmark.peaks(peakSet="narrow",
                                peaks=atac.allpeaks.pbmc,
                                peakcalls=(peakcalls.pbmc %>% 
                                             mutate(peakco=paste(chr,start,end,sep = "_")) %>%
                                             data.frame(.,row.names = "peakco",stringsAsFactors = F))[atac.allpeaks.pbmc$peakco,],
                                anno=atac.allpeaks.anno.pbmc,
                                samps=sampid,
                                libsize=libsize.atac.pbmc,
                                depth=depth,
                                group=age.group,
                                max.logFC.quantile=0.25,
                                min.maxct.quantile=0.9,
                                min.peakcount.quantile=0.9,
                                filter.thresh=0.925,
                                show.plots=TRUE)