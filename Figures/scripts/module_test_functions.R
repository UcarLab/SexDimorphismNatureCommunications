# Helper functions for aging paper enrichment analysis

# Module enrichment
module_test <- function(peak.anno=NULL,
                        background=NULL,
                        test.fdr=0.05,
                        test.logfc=0,
                        logfc.name="logFC",
                        fdr.name="FDR",
                        direction.string=list(pos="Opening with age",neg="Closing with age"),
                        comp.name=NULL,
                        geneset.path=NULL,
                        enrich.fdr=0.05,
                        flank.size=5e4,
                        geneset.sel="vp2008",
                        synonymize.ext=FALSE,
                        verbose=FALSE) {
  
  # require(GenomicRanges)
  require(dplyr)
  
  genesets <- load(geneset.path)
  geneset.sel <- geneset.sel[1]
  rm(list = genesets[grep(geneset.sel,genesets,invert = T)])
  geneset.genes <- get(paste("geneset.genes.",geneset.sel,sep = ""))
  geneset.names <- get(paste("geneset.names.",geneset.sel,sep = ""))
  geneset.name <- get(paste("geneset.name.",geneset.sel,sep = ""))
  min.modcount <- get(paste("min.modcount.",geneset.sel,sep = ""))
  min.genecount <- get(paste("min.genecount.",geneset.sel,sep = ""))

  if (is.null(background)) {
    print("Warning: no background specified. Using the list of genes in the input set.")
    background <- sort(unique(as.character(peak.anno$GeneName)))
  }
  peak.anno_filtered <- peak.anno %>%
    rename_(logFC=logfc.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr|abs(logFC) >= test.logfc) %>%
    mutate(direction=ifelse(logFC<0,direction.string$neg,direction.string$pos)) %>%
    filter(abs(DistancetoTSS)<=flank.size)
  
  if (nrow(peak.anno_filtered)==0) return(NULL)
  
  # ss <- GRanges(peak.anno_filtered)
  
  # if (synonymize.ext) {
  #   ss@elementMetadata$GeneName <- sub("-IT.$","",sub("-AS.$","",ss@elementMetadata$GeneName))
  #   geneset.genes$GeneName <- sub("-IT.$","",sub("-AS.$","",geneset.genes$GeneName))
  # }
  # names(ss@elementMetadata) <- sub(paste(".",comp.name,sep=""),"",names(ss@elementMetadata))
  # browser()
  geneset.genes.do <- geneset.genes[is.element(geneset.genes$GeneName,background),]
  geneset.do <- list(genes=geneset.genes.do,names=geneset.names,labels=geneset.name)
  
  # dir.test <- split(ss,f = ss@elementMetadata$direction)
  dir.test <- split(peak.anno_filtered,f = peak.anno_filtered$direction)
  if (!any(unlist(lapply(dir.test,function(X) length(X)>=min.genecount)))) return(NULL)
  
  genesets.test <- lapply(dir.test[unlist(lapply(dir.test,length))>0], function(X) {
    # gL = unique(as.character(X@elementMetadata$GeneName))
    gL = unique(as.character(X$GeneName))
    if (length(gL)!=0) sort(gL[gL %in% background])
  })
  genesets.test <- genesets.test[unlist(lapply(genesets.test,function(X) length(X)>=min.genecount))]
  
  gsets.do <- list(genesets=genesets.test,labels=genesets.test,background=background)
  moduletest.test <- modules.enrichtest(geneset.do,min.modcount,min.genecount,enrich.fdr,gsets.do)
  if (is.null(moduletest.test)) return(NULL)
  moduletest.test$hypergeom.fdr <- p.adjust(10^-(moduletest.test$hypergeom.p-1e-17),method = "fdr")
  
  moduletest.test.top <- moduletest.test[moduletest.test$Module.Name %in% unique(as.character(moduletest.test[moduletest.test$hypergeom.fdr<enrich.fdr,]$Module.Name)),]
  moduletest.out <- list(all=moduletest.test,top=moduletest.test.top,geneset=geneset.do)
  
  if (verbose) print("Done.")
  
  # detach("package:GenomicRanges", unload=TRUE)
  return(moduletest.out)
}

module_test_rna <- function(rna.glm=NULL,
                            background=NULL,
                            test.fdr=0.05,
                            test.logfc=0,
                            logfc.name="logFC",
                            fdr.name="FDR",
                            direction.string=list(pos="Opening with age",neg="Closing with age"),
                            comp.name=NULL,
                            geneset.path=NULL,
                            enrich.fdr=0.05,
                            geneset.sel="vp2008",
                            synonymize.ext=FALSE,
                            verbose=FALSE) {
  
  require(dplyr)
  
  genesets <- load(geneset.path)
  geneset.sel <- geneset.sel[1]
  rm(list = genesets[grep(geneset.sel,genesets,invert = T)])
  geneset.genes <- get(paste("geneset.genes.",geneset.sel,sep = ""))
  geneset.names <- get(paste("geneset.names.",geneset.sel,sep = ""))
  geneset.name <- get(paste("geneset.name.",geneset.sel,sep = ""))
  min.modcount <- get(paste("min.modcount.",geneset.sel,sep = ""))
  min.genecount <- get(paste("min.genecount.",geneset.sel,sep = ""))
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of genes in the input set.")
    background <- sort(unique(as.character(rna.glm$GeneName)))
  }
  rna.glm_filtered <- rna.glm %>%
    rename_(logFC=logfc.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr|abs(logFC) >= test.logfc) %>%
    mutate(direction=ifelse(logFC<0,direction.string$neg,direction.string$pos))
  
  if (nrow(rna.glm_filtered)==0) return(NULL)
  
  if (synonymize.ext) {
    rna.glm_filtered$GeneName <- sub("-IT.$","",sub("-AS.$","",rna.glm_filtered$GeneName))
    geneset.genes$GeneName <- sub("-IT.$","",sub("-AS.$","",geneset.genes$GeneName))
  }
  
  geneset.genes.do <- geneset.genes[is.element(geneset.genes$GeneName,background),]
  geneset.do <- list(genes=geneset.genes.do,names=geneset.names,labels=geneset.name)
  
  dir.test <- split(rna.glm_filtered,f = rna.glm_filtered$direction)
  if (!any(unlist(lapply(dir.test,function(X) length(X)>=min.genecount)))) return(NULL)
  
  genesets.test <- lapply(dir.test[unlist(lapply(dir.test,length))>0], function(X) {
    gL = unique(as.character(X$GeneName))
    if (length(gL)!=0) sort(gL[gL %in% background])
  })
  genesets.test <- genesets.test[unlist(lapply(genesets.test,function(X) length(X)>=min.genecount))]
  
  gsets.do <- list(genesets=genesets.test,labels=genesets.test,background=background)
  moduletest.test <- modules.enrichtest(geneset.do,min.modcount,min.genecount,enrich.fdr,gsets.do)
  if (is.null(moduletest.test)) return(NULL)
  moduletest.test$hypergeom.fdr <- p.adjust(10^-(moduletest.test$hypergeom.p-1e-17),method = "fdr")
  
  moduletest.test.top <- moduletest.test[moduletest.test$Module.Name %in% unique(as.character(moduletest.test[moduletest.test$hypergeom.fdr<enrich.fdr,]$Module.Name)),]
  moduletest.out <- list(all=moduletest.test,top=moduletest.test.top,geneset=geneset.do)
  
  if (verbose) print("Done.")
  
  return(moduletest.out)
}

module_test_peaks <- function(peak.anno=NULL,
                              background=NULL,
                              test.fdr=0.05,
                              test.logfc=0,
                              logfc.name="logFC",
                              fdr.name="FDR",
                              direction.string=list(pos="Opening with age",neg="Closing with age"),
                              comp.name=NULL,
                              peakset=NULL,
                              min.peakcount=1,
                              enrich.fdr=0.05,
                              flank.size=5e4,
                              verbose=FALSE) {
  
  require(GenomicRanges)
  require(dplyr)
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of peaks in the input set.")
    background <- peak.anno %>% select(chr,start,end) %>% unique()
  }
  peak.anno_filtered <- peak.anno %>%
    rename_(logFC=logfc.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr|abs(logFC) >= test.logfc) %>%
    mutate(direction=ifelse(logFC<0,direction.string$neg,direction.string$pos)) %>%
    filter(abs(DistancetoTSS)<=flank.size)
  
  if (nrow(peak.anno_filtered)==0) return(NULL)
  
  peaks <- GRanges(peak.anno_filtered)
  bg <- GRanges(background)
  peaksets <- GRangesList(lapply(peakset,GRanges))
  
  dir.test <- split(peaks,f=peaks$direction)
  if (!any(unlist(lapply(dir.test,function(X) length(X)>=min.peakcount)))) return(NULL)
  
  chstate.test <- do.call(rbind,lapply(names(dir.test), function(n) {
    pks <- dir.test[[n]]
    do.call(rbind,lapply(names(peaksets), function(ct) {
      data.frame(celltype=ct,
                 hypergeom.p=-log10(phyper(length(intersect(pks,peaksets[[ct]])),length(intersect(bg,peaksets[[ct]])),length(setdiff(bg,peaksets[[ct]])),length(pks),lower.tail = F)),
                 peak.count=length(intersect(pks,peaksets[[ct]])),
                 direction=n,
                 stringsAsFactors = F)
    }))
  })) %>%
    mutate(signed.hypergeom.p=ifelse(grepl("Closing",direction),-hypergeom.p,hypergeom.p),
           hypergeom.fdr=p.adjust(10^(-hypergeom.p),method = "fdr"))
  
  chstate.test.top <- chstate.test[chstate.test$celltype %in% unique(as.character(chstate.test[chstate.test$hypergeom.fdr<enrich.fdr,]$celltype)),]
  chstatetest.out <- list(all=chstate.test,top=chstate.test.top,peaksets=peaksets)
  
  if (verbose) print("Done.")
  
  # detach("package:GenomicRanges", unload=TRUE)
  return(chstatetest.out)
}

atac.module_test_clusters <- function(peak.anno=NULL,
                                 background=NULL,
                                 test.fdr=0.05,
                                 clusters.name="Cluster",
                                 fdr.name="FDR",
                                 comp.name=NULL,
                                 geneset.path=NULL,
                                 enrich.fdr=0.05,
                                 flank.size=5e4,
                                 geneset.sel=c("vp2008","wp"),
                                 synonymize.ext=FALSE,
                                 verbose=FALSE) {
  
  require(GenomicRanges)
  require(dplyr)
  
  geneset.sel <- geneset.sel[1]
  
  genesets <- load(geneset.path)
  rm(list = genesets[grep(geneset.sel,genesets,invert = T)])
  geneset.genes <- get(paste("geneset.genes.",geneset.sel,sep = ""))
  geneset.names <- get(paste("geneset.names.",geneset.sel,sep = ""))
  geneset.name <- get(paste("geneset.name.",geneset.sel,sep = ""))
  min.modcount <- get(paste("min.modcount.",geneset.sel,sep = ""))
  min.genecount <- get(paste("min.genecount.",geneset.sel,sep = ""))
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of genes in the input set.")
    background <- sort(unique(as.character(peak.anno$GeneName)))
  }
  peak.anno_filtered <- peak.anno %>%
    rename_(Cluster=clusters.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr) %>%
    filter(abs(DistancetoTSS)<=flank.size)
  
  if (nrow(peak.anno_filtered)==0) return(NULL)
  
  ss <- GRanges(peak.anno_filtered)
  
  if (synonymize.ext) {
    ss@elementMetadata$GeneName <- sub("-IT.$","",sub("-AS.$","",ss@elementMetadata$GeneName))
    geneset.genes$GeneName <- sub("-IT.$","",sub("-AS.$","",geneset.genes$GeneName))
  }
  names(ss@elementMetadata) <- sub(paste(".",comp.name,sep=""),"",names(ss@elementMetadata))
  
  geneset.genes.do <- geneset.genes[is.element(geneset.genes$GeneName,background),]
  geneset.do <- list(genes=geneset.genes.do,names=geneset.names,labels=geneset.name)
  
  clust.test <- split(ss,f = ss@elementMetadata$Cluster)
  if (!any(unlist(lapply(clust.test,function(X) length(X)>=min.genecount)))) return(NULL)
  
  genesets.test <- lapply(clust.test[unlist(lapply(clust.test,length))>0], function(X) {
    gL = unique(as.character(X@elementMetadata$GeneName))
    if (length(gL)!=0) sort(gL[gL %in% background])
  })
  genesets.test <- genesets.test[unlist(lapply(genesets.test,function(X) length(X)>=min.genecount))]
  
  gsets.do <- list(genesets=genesets.test,labels=genesets.test,background=background)
  moduletest.test <- modules.enrichtest(geneset.do,min.modcount,min.genecount,enrich.fdr,gsets.do)
  if (is.null(moduletest.test)) return(NULL)
  moduletest.test$hypergeom.fdr <- p.adjust(10^-(moduletest.test$hypergeom.p-1e-17),method = "fdr")
  
  moduletest.test.top <- moduletest.test[moduletest.test$Module.Name %in% unique(as.character(moduletest.test[moduletest.test$hypergeom.fdr<enrich.fdr,]$Module.Name)),]
  moduletest.out <- list(all=moduletest.test,top=moduletest.test.top,geneset=geneset.do)
  
  if (verbose) print("Done.")
  
  return(moduletest.out)
}

atac.chrstates_test_clusters <- function(peak.anno=NULL,
                                         background=NULL,
                                         clusters.name="Cluster",
                                         comp.name=NULL,
                                         peakset=NULL,
                                         min.peakcount=1,
                                         enrich.fdr=0.05,
                                         flank.size=5e4,
                                         verbose=FALSE) {
  
  # require(GenomicRanges)
  require(dplyr)
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of peaks in the input set.")
    background <- peak.anno %>% select(chr,start,end) %>% unique()
  }
  peak.anno_filtered <- peak.anno %>%
    filter(abs(DistancetoTSS)<=flank.size) %>%
    rename(Cluster=clusters.name)
  
  if (nrow(peak.anno_filtered)==0) return(NULL)
  
  peaks <- GenomicRanges::GRanges(peak.anno_filtered)
  bg <- GenomicRanges::GRanges(background)
  peaksets <- GenomicRanges::GRangesList(lapply(peakset,GenomicRanges::GRanges))
  
  dir.test <- split(peaks,f=peaks$Cluster)
  if (!any(unlist(lapply(dir.test,function(X) length(X)>=min.peakcount)))) return(NULL)
  
  chstate.test <- do.call(rbind,lapply(names(dir.test), function(n) {
    pks <- dir.test[[n]]
    do.call(rbind,lapply(names(peaksets), function(ct) {
      data.frame(celltype=ct,
                 hypergeom.p=-log10(phyper(length(intersect(pks,peaksets[[ct]])),length(intersect(bg,peaksets[[ct]])),length(setdiff(bg,peaksets[[ct]])),length(pks),lower.tail = F)),
                 peak.count=length(intersect(pks,peaksets[[ct]])),
                 Cluster=n,
                 stringsAsFactors = F)
    }))
  })) %>%
    mutate(hypergeom.fdr=p.adjust(10^(-hypergeom.p),method = "fdr"))
  
  chstate.test.top <- chstate.test[chstate.test$celltype %in% unique(as.character(chstate.test[chstate.test$hypergeom.fdr<enrich.fdr,]$celltype)),]
  chstatetest.out <- list(all=chstate.test,top=chstate.test.top,peaksets=peaksets)
  
  if (verbose) print("Done.")
  
  # detach("package:GenomicRanges", unload=TRUE)
  return(chstatetest.out)
}

atac.chrstates_module_test <- function(peak.anno=NULL,
                                       background=NULL,
                                       test.fdr=0.05,
                                       test.logfc=0,
                                       logfc.name="logFC",
                                       fdr.name="FDR",
                                       direction.string=list(pos="Opening with age",neg="Closing with age"),
                                       comp.name=NULL,
                                       peakset=NULL,
                                       min.peakcount=1,
                                       enrich.fdr=0.05,
                                       flank.size=5e4,
                                       verbose=FALSE) {
  
  require(GenomicRanges)
  require(dplyr)
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of peaks in the input set.")
    background <- peak.anno %>% select(chr,start,end) %>% unique()
  }
  peak.anno_filtered <- peak.anno %>%
    rename_(logFC=logfc.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr|abs(logFC) >= test.logfc) %>%
    mutate(direction=ifelse(logFC<0,direction.string$neg,direction.string$pos)) %>%
    filter(abs(DistancetoTSS)<=flank.size)
  
  if (nrow(peak.anno_filtered)==0) return(NULL)
  
  peaks <- GRanges(peak.anno_filtered)
  bg <- GRanges(background)
  peaksets <- GRangesList(lapply(peakset,GRanges))
  
  dir.test <- split(peaks,f=peaks$direction)
  if (!any(unlist(lapply(dir.test,function(X) length(X)>=min.peakcount)))) return(NULL)
  
  chstate.test <- do.call(rbind,lapply(names(dir.test), function(n) {
    pks <- dir.test[[n]]
    do.call(rbind,lapply(names(peaksets), function(ct) {
      data.frame(celltype=ct,
                 hypergeom.p=-log10(phyper(length(intersect(pks,peaksets[[ct]])),length(intersect(bg,peaksets[[ct]])),length(setdiff(bg,peaksets[[ct]])),length(pks),lower.tail = F)),
                 peak.count=length(intersect(pks,peaksets[[ct]])),
                 direction=n,
                 stringsAsFactors = F)
    }))
  })) %>%
    mutate(signed.hypergeom.p=ifelse(grepl("Closing",direction),-hypergeom.p,hypergeom.p),
           hypergeom.fdr=p.adjust(10^(-hypergeom.p),method = "fdr"))
  
  chstate.test.top <- chstate.test[chstate.test$celltype %in% unique(as.character(chstate.test[chstate.test$hypergeom.fdr<enrich.fdr,]$celltype)),]
  chstatetest.out <- list(all=chstate.test,top=chstate.test.top,peaksets=peaksets)
  
  if (verbose) print("Done.")
  
  # detach("package:GenomicRanges", unload=TRUE)
  return(chstatetest.out)
}

rna.module_test_clusters <- function(rna.anno=NULL,
                                     background=NULL,
                                     test.fdr=0.05,
                                     clusters.name="Cluster",
                                     fdr.name="FDR",
                                     comp.name=NULL,
                                     geneset.path=NULL,
                                     enrich.fdr=0.05,
                                     geneset.sel=c("vp2008","wp"),
                                     synonymize.ext=FALSE,
                                     verbose=FALSE) {
  
  require(dplyr)
  
  geneset.sel <- geneset.sel[1]
  
  genesets <- load(geneset.path)
  rm(list = genesets[grep(geneset.sel,genesets,invert = T)])
  geneset.genes <- get(paste("geneset.genes.",geneset.sel,sep = ""))
  geneset.names <- get(paste("geneset.names.",geneset.sel,sep = ""))
  geneset.name <- get(paste("geneset.name.",geneset.sel,sep = ""))
  min.modcount <- get(paste("min.modcount.",geneset.sel,sep = ""))
  min.genecount <- get(paste("min.genecount.",geneset.sel,sep = ""))
  
  if (is.null(background)) {
    print("Warning: no background specified. Using the list of genes in the input set.")
    background <- sort(unique(as.character(rna.anno.clust$GeneName)))
  }
  rna.anno_filtered <- rna.anno %>%
    rename_(Cluster=clusters.name,FDR=fdr.name) %>%
    filter(FDR <= test.fdr)
  
  if (nrow(rna.anno_filtered)==0) return(NULL)
  
  if (synonymize.ext) {
    rna.anno_filtered$GeneName <- sub("-IT.$","",sub("-AS.$","",rna.anno_filtered$GeneName))
    geneset.genes$GeneName <- sub("-IT.$","",sub("-AS.$","",geneset.genes$GeneName))
  }
  
  geneset.genes.do <- geneset.genes[is.element(geneset.genes$GeneName,background),]
  geneset.do <- list(genes=geneset.genes.do,names=geneset.names,labels=geneset.name)
  
  clust.test <- split(rna.anno_filtered,f = rna.anno_filtered$Cluster)
  if (!any(unlist(lapply(clust.test,function(X) length(X)>=min.genecount)))) return(NULL)
  
  genesets.test <- lapply(clust.test[unlist(lapply(clust.test,length))>0], function(X) {
    gL = unique(as.character(X$GeneName))
    if (length(gL)!=0) sort(gL[gL %in% background])
  })
  genesets.test <- genesets.test[unlist(lapply(genesets.test,function(X) length(X)>=min.genecount))]
  
  gsets.do <- list(genesets=genesets.test,labels=genesets.test,background=background)
  moduletest.test <- modules.enrichtest(geneset.do,min.modcount,min.genecount,enrich.fdr,gsets.do)
  moduletest.test$hypergeom.fdr <- p.adjust(10^-(moduletest.test$hypergeom.p-1e-17),method = "fdr")
  
  moduletest.test.top <- moduletest.test[moduletest.test$Module.Name %in% unique(as.character(moduletest.test[moduletest.test$hypergeom.fdr<enrich.fdr,]$Module.Name)),]
  moduletest.out <- list(all=moduletest.test,top=moduletest.test.top,geneset=geneset.do)
  
  if (verbose) print("Done.")
  
  return(moduletest.out)
}
