# Sex x aging study
# Module: modules.enrichtest.R
# Description: Function: accepts lists of genes and a module family, and calculates enrichment of each module or pathway or gene set in that family

library(ggplot2)
library(reshape2)
library(MKmisc)
library(NMF)

sse <- function(x) sum((x-mean(x))^2)
andred <- function(x) {Reduce("&",x)}
orred <- function(x) {Reduce("|",x)}

modules.enrichtest <- function(modset,min.modcount,min.genecount,fdr.cut,gsets) {
  
  modset.name <- modset[["labels"]]
  modset.genes <- modset[["genes"]]
  modset.names <- modset[["names"]]
  
  genesets <- gsets[["genesets"]]
  genesets.labels <- gsets[["labels"]]
  genesets.sign <- gsets[["sign"]]
  gene.uni <- gsets[["background"]]
  
  modset.modules <- merge(modset.genes,modset.names, by = "Module.ID")
  module.counts <- data.frame(table(modset.modules$Module.Name))
  colnames(module.counts) <- c("Module.Name","Module.Counts")
  module.counts$Module.in <- module.counts$Module.Counts>=min.modcount
  
  # Update modules
  modset.modules <- merge(modset.modules,module.counts,by = "Module.Name")
  modset.modules <- modset.modules[modset.modules$Module.in,]
  modset.modules <- modset.modules[,-which(is.element(colnames(modset.modules), c("Module.Counts","Module.in")))]
  modset.modules$Module.Name <- factor(modset.modules$Module.Name)
  modset.modules.genesets <- unique(modset.modules[,c("Module.Name","GeneName")])
  gene.mod <- as.vector(unique(modset.modules.genesets$GeneName)) # genes annotated to a module/pathway
  
  # "enrichment" table
  modset.counts <- data.frame(Total=table(modset.modules.genesets$Module.Name), row.names = 1)
  modset.counts.text <- data.frame(row.names=rownames(modset.counts))
  modset.phyper <- data.frame(row.names = rownames(modset.counts))
  modset.prop <- data.frame(row.names = rownames(modset.counts))
  modset.jaccard <- data.frame(row.names = rownames(modset.counts))
  modset.cellnote <- data.frame(row.names = rownames(modset.counts))
  for (i in 1:length(genesets)) {
    gene.sampler <- intersect(genesets[[i]],gene.uni)
    modset.counts <- data.frame(merge(modset.counts,data.frame(setNames(data.frame(table(merge(modset.modules.genesets,data.frame(GeneName=genesets[[i]]))$Module.Name)),c("Total.Freq",names(genesets)[i])),check.names = F,row.names = 1),by = "row.names"),row.names = 1,check.names = F)
    modset.counts.text <- data.frame(merge(modset.counts.text,setNames(data.frame(new=apply(modset.counts[,c(names(genesets)[i],"Total.Freq")],1,function(x) paste(x[1],"/",x[2], sep = ""))),names(genesets)[i]), by= "row.names"), row.names = 1,check.names = F)
    modset.phyper <- data.frame(merge(modset.phyper,setNames(data.frame(new=apply(modset.counts[,c(names(genesets)[i],"Total.Freq")],1,function(x) -log10(phyper(x[1], x[2], length(gene.uni)-x[2], length(gene.sampler), lower.tail = F)))),names(genesets)[i]), by= "row.names"), row.names = 1,check.names = F)
    modset.prop <- data.frame(merge(modset.prop,setNames(data.frame(new=apply(modset.counts[,c(names(genesets)[i],"Total.Freq")],1,function(x) x[1]/x[2])),names(genesets)[i]), by= "row.names"), row.names = 1,check.names = F)
    modset.jaccard <- data.frame(merge(modset.jaccard,setNames(data.frame(new=apply(modset.counts[,c(names(genesets)[i],"Total.Freq")],1,function(x) x[1]/(x[2]+length(genesets[[i]])-x[1]))),names(genesets)[i]), by= "row.names"), row.names = 1,check.names = F)
    modset.cellnote <- data.frame(merge(modset.cellnote,setNames(data.frame(new=apply(modset.counts[,c(names(genesets)[i],"Total.Freq")],1,
                                                                                      function(x) paste("Count=",x[1],"/",x[2],
                                                                                                        " Jaccard=",sprintf("%g",x[1]/(x[2]+length(genesets[[i]])-x[1])),
                                                                                                        " loginvP.hyper=",sprintf("%g",-log10(phyper(x[1], x[2], length(gene.uni)-x[2], length(gene.sampler), lower.tail = F))),
                                                                                                        sep = "")
    )
    ),names(genesets)[i]),
    by= "row.names"),
    row.names = 1,check.names = F)
  }
  
  # clean up table
  minlogP <- -log10(fdr.cut) # remove modules/pahways for which no gene set reaches this loginvP
  minCount <- min.genecount # remove modules/pathways for which no gene set reaches this count
  minSelect <- rowSums(setNames(data.frame(modset.counts[,-1]),names(modset.counts)[-1])>=minCount)>0 | rowSums(modset.phyper>minlogP)>0
  
  if (!any(minSelect)) return(NULL)
  
  modset.phyper <- subset(modset.phyper,minSelect)
  modset.counts.text <- subset(modset.counts.text,minSelect)
  modset.counts <- subset(modset.counts,minSelect)
  modset.prop <- subset(modset.prop,minSelect)
  modset.jaccard <- subset(modset.jaccard,minSelect)
  modset.cellnote <- subset(modset.cellnote,minSelect)
  
  modset.phyper.glyph <- modset.phyper
  modset.phyper.glyph$Module.Name <- rownames(modset.phyper.glyph)
  rownames(modset.phyper.glyph) <- NULL
  modset.phyper.glypht <- melt(modset.phyper.glyph,id.vars = "Module.Name", variable.name = "gene.set", value.name = "hypergeom.p")
  modset.phyper.glypht$hypergeom.p <- as.vector(modset.phyper.glypht$hypergeom.p)
  
  modset.counts.glyph <- modset.counts
  modset.counts.glyph$Module.Name <- rownames(modset.counts.glyph)
  rownames(modset.counts.glyph) <- NULL
  modset.counts.preglyph <- modset.counts.glyph
  modset.counts.preglyph[,-c(1,ncol(modset.counts.preglyph))] <- modset.counts.preglyph[,1]-modset.counts.preglyph[,-c(1,ncol(modset.counts.preglyph))]
  modset.counts.glypht <- rbind(cbind(melt(modset.counts.glyph[,-1],id.vars = "Module.Name", variable.name = "gene.set", value.name = "gene.count"),dummypos=1),
                                cbind(melt(modset.counts.preglyph[,-1],id.vars = "Module.Name", variable.name = "gene.set", value.name = "gene.count"),dummypos=2))
  modset.counts.glypht$gene.count <- as.vector(modset.counts.glypht$gene.count)
  modset.counts.glypht$gene.set <- as.vector(modset.counts.glypht$gene.set)
  modset.counts.glypht <- modset.counts.glypht[order(modset.counts.glypht$Module.Name,modset.counts.glypht$gene.set,modset.counts.glypht$dummypos),]
  modset.counts.glypht <- merge(merge(modset.counts.glypht,modset.counts.preglyph[,c("Total.Freq","Module.Name")]),modset.phyper.glypht)
  modset.counts.glypht$Module.Label <- paste(modset.counts.glypht$Module.Name," (",modset.counts.glypht$Total.Freq," genes)", sep = "")
  
  return(modset.counts.glypht[modset.counts.glypht$dummypos==1,])
}