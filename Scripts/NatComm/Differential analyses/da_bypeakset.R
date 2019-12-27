# Function to run differential analysis (DA) on ATAC-seq for the sex:age data using the specified models and parameters.
# Function assumes ATACseq data objects obtained after benchmarking are loaded in the working environment.
# Example call used in these analyses:
# 
# # Global parameters used in DA
# global.pars <- list(use.adjustment="libsize", # options: no | libsize | frips | depth
#                     use.combat=T,
#                     batch.variable="batchdate",
#                     use.sva=F,
#                     n.sva="auto",
#                     use.robust=T,
#                     use.tagwise=T,
#                     norm.method="TMM",
#                     trend.method="locfit",
#                     min.peakcalls=3,
#                     max.dist2tss=1e50,
#                     run.name=run.name)
# if (!global.pars$use.sva) global.pars$n.sva=0
# 
# design <- paste("model.matrix(~sex +",global.pars$batch.variable,"+ sex:age.group)")
# design.null <- paste0("model.matrix(~",global.pars$batch.variable,")")
# 
# # Data frame with relevant factors, individually loaded from filtered_global.stats_narrowPeaks_filtered.RData
# metadata.factors <- data.frame(sex=as.character(sex),
#                                age=age,
#                                age.group=as.character(age.group),
#                                batch=batch,
#                                batchdate=as.character(batchdate),
#                                season=season,
#                                sampid=as.character(samps),
#                                nreps=nreps,
#                                depth=log2(depth),
#                                libsize=log2(libsize.atac.pbmc),
#                                frips=frips.atac.pbmc,
#                                stringsAsFactors = F)
# 
# # Basic model is adjusted depending on whether library size or other factors is included as covariate, and
# # whether ComBat is used to correct for batch effects prior to running the analyses
# if (global.pars$use.combat) {
#   if (global.pars$use.adjustment!="no") {
#     design <- gsub(global.pars$batch.variable,global.pars$use.adjustment,design)
#     design.null <- gsub(global.pars$batch.variable,global.pars$use.adjustment,design.null)
#   } else {
#     design <- gsub(paste(global.pars$batch.variable,"\\+ "),"",design)
#     design.null <- gsub(paste(global.pars$batch.variable,"\\+ "),"",design.null)
#   }
# } else {
#   if (global.pars$use.adjustment!="no") {
#     design <- gsub(global.pars$batch.variable,paste(global.pars$batch.variable,"+",global.pars$use.adjustment),design)
#     design.null <- gsub(global.pars$batch.variable,paste(global.pars$batch.variable,"+",global.pars$use.adjustment),design.null)
#   }
# }
# atac.design <- with(metadata.factors,eval(parse(text=design)))
# atac.design.null <- with(metadata.factors,eval(parse(text=design.null)))
# 
# # Batch correction via ComBat
# peaks <- atac.allpeaks.pbmc
# if (globalPars$use.combat) {
#   cat("Applying batch correction via ComBat, batch variable:",globalPars$batch.variable)
#   library(sva)
#   peaks.combadj <- expm1(ComBat(dat = as.matrix(log1p(peaks[,metadata.factors$sampid])),batch = get(globalPars$batch.variable)[metadata.factors$sampid],mod = atac.design))
#   peaks.combadj[peaks.combadj < 0] <- 0 # zero negative values
#   peaks.combadj <- setNames(data.frame(as.matrix(peaks.combadj) %*% diag(libsize/colSums(peaks.combadj))),names(libsize)) # scale to original library sizes
#   peaks <- cbind(peaks %>% select(chr,start,end,peakco,nhits),peaks.combadj)
# }
# 
# # Deinition of contrasts to test, variable whether DA compare age groups by sex, or sexes by age group
# # This is for analysis by sex
# focus.contrast.factors = "age"
# focus.contrast.set <- c(as.list(colnames(atac.design)[grep(paste(focus.contrast.factors,collapse="|"),colnames(atac.design))]),
#                         list(c(as.numeric(grepl("sexF\\:",colnames(atac.design))),rep(0,ncol(da.list$design)-ncol(atac.design))),
#                              c(as.numeric(grepl("sexM\\:",colnames(atac.design))),rep(0,ncol(da.list$design)-ncol(atac.design)))))
# focus.contrast.set <- lapply(focus.contrast.set,function(x) {
#   if (is.numeric(x)) {
#     x[min(which(as.logical(x)))] <- -1
#   }
#   return(x)
# })
# focus.contrast.set.names <- sapply(focus.contrast.set, function(x) {
#   if (is.numeric(x)) {
#     n <- sub("sex.*age","Age",sub("^sexM\\:age.groupage","Males_Age",sub("^sexF\\:age.groupage","Females_Age",paste(rev(colnames(atac.design)[as.logical(x)]),collapse="x"))))
#   } else {
#     n <- paste0(sub("M","Males_Age",sub("F","Females_Age",sub("sex","",sub("\\_.*age","",sub("\\:","_",x))))),"xAge1")
#   }
# })
# 
# # Differential analysis routine call
# da_bysex <- lapply(setNames(focus.contrast.set,focus.contrast.set.names), function(focus.contrast) {
#   cat("Now processing:",focus.contrast,"...\n")
#   da.list <- da_bypeakset(atac.data.raw=peaks,
#                           atac.peakcalls=peakcalls.pbmc,
#                           atac.anno=atac.allpeaks.anno.pbmc,
#                           atac.design=atac.design,
#                           atac.design.null=atac.design.null,
#                           focus.contrast=focus.contrast,
#                           samp.names=samps,
#                           run.name=focus.contrast,
#                           min.peakcalls=3,
#                           max.dist2tss=5e10,
#                           group.names=age.group,
#                           atac.peakset=toupper(peak.set),
#                           atac.libsize=libsize.atac.pbmc,
#                           use.adjustment="no",
#                           use.sva=T,
#                           n.sva="auto",
#                           use.robust=T,
#                           use.tagwise=T,
#                           norm.method="TMM",
#                           trend.method="locfit")
#   cat("Found",sum(da.list$glm$FDR<=0.05),"hits at",100*0.05,"% FDR\n")
#   da.output <- list(da.list=da.list,
#                     metadata.factors=metadata.factors,
#                     atac.design=atac.design,
#                     atac.design.null=atac.design.null,
#                     focus.contrast=focus.contrast,
#                     global.pars=global.pars,
#                     age.group=age.group)
#   return(da.output)
# })
#################################################################################################################################################################


da_bypeakset <- function(atac.dge=NULL, # if desired, a DGE object from an earlier run of EdgeR can be used to retrieve data and analysis parameters
                         atac.data.raw=NULL, # ignored if atac.dge present
                         atac.peakcalls=NULL, # ignored if atac.dge present
                         atac.anno=NULL, # ignored if atac.dge present
                         atac.peakset=NULL, # ignored if atac.dge present
                         atac.design=NULL, # ignored if atac.dge present
                         atac.design.null=NULL,
                         focus.contrast=NULL,
                         atac.libsize=NULL,
                         samp.names=NULL, # ignored if atac.dge present
                         group.names=NULL, # ignored if atac.dge present
                         run.name=NULL, # ignored if atac.dge present
                         min.peakcalls=1, # ignored if atac.dge present
                         max.dist2tss=5e4, # ignored if atac.dge present
                         exclude.samples=NULL, # ignored if atac.dge present
                         use.adjustment=c("no","libsize"),
                         use.sva=F, # ignored if atac.dge present
                         n.sva="auto", # ignored if atac.dge present
                         use.robust=T, # ignored if atac.dge present
                         use.tagwise=F, # ignored if atac.dge present
                         norm.method=c("TMM","RLE","upperquartile","cpm","none"), # ignored if atac.dge present
                         trend.method="none") { # ignored if atac.dge present
  require(edgeR)
  require(sva)
  require(reshape2)
  require(dplyr)
  
  if (is.null(atac.dge)) {
    if (!is.null(exclude.samples)) {
      samp.names <- setdiff(samp.names,exclude.samples)
    }
    atac.libsize <- atac.libsize[samp.names]
    group.names <- droplevels(group.names[samp.names])
    
    n_init <- nrow(atac.data.raw)
    atac.data.raw <- atac.data.raw %>% 
      select(c("chr","start","end","nhits","peakco",samp.names)) %>%
      left_join(atac.anno %>% select(peakco,DistancetoTSS),by = "peakco") %>%
      left_join(atac.peakcalls[,c("chr","start","end",samp.names)] %>% 
                  slice(match(rownames(.),rownames(atac.data.raw))) %>%
                  mutate(Total.counts=rowSums(.[,samp.names])) %>% 
                  select(chr,start,end,Total.counts) %>% 
                  data.frame(.,row.names = rownames(atac.data.raw)),
                by = c("chr","start","end")) %>%
      filter(!abs(DistancetoTSS)>max.dist2tss,!Total.counts<min.peakcalls)
    
    atac.mx.raw <- as.matrix(atac.data.raw[,samp.names])
    rownames(atac.mx.raw) <- atac.data.raw$peakco
    
    cat(paste("Run name:",run.name,"\n"))
    cat(paste("Input data:",deparse(quote(atac.data.raw)),"\n"))
    cat(paste("Input peak set:",atac.peakset,"\n"))
    cat(paste("Sample names:",paste(samp.names,collapse = "-"),"\n"))
    cat(paste("group names:",paste(levels(group.names),collapse = ","),"\n"))
    cat(paste("Peak filter 1: each peak must be called at least",min.peakcalls,"time(s) to enter the analysis.\n"))
    cat(paste("Peak filter 2: gene annotations must lie at a distance of at most",max.dist2tss,"bp to enter the analysis.\n"))
    cat(paste("==>",nrow(atac.data.raw),"out of",n_init,"peaks entered the analysis.\n"))
    cat(paste("Adjusted: ",use.adjustment,"\n"))
    cat(paste("Normalization method: ",norm.method,"\n"))
    cat(paste("Using Tagwise dispersion estimates:",use.tagwise,"\n"))
    cat(paste("Dispersion estimation trend method:",trend.method,"\n"))
    cat(paste("Using SVA:",use.sva,"\n"))
    if (use.sva) cat(paste("N SV:",n.sva,"\n"))
    cat(paste("Using robust estimate of dispersion factors:",use.robust,"\n"))
    
    cat("Building data matrix...\n")
    if (norm.method!="cpm") {
      atac.dge <- DGEList(atac.mx.raw,group=group.names,lib.size = atac.libsize)
      cat(paste("Applying",norm.method,"normalization...\n"))
      atac.dge <- calcNormFactors(atac.dge,method = norm.method);
      atac.mx.norm <- atac.mx.raw %*% diag(atac.dge$samples$norm.factors)
    } else {
      cat(paste("Applying",norm.method,"normalization...\n"))
      atac.mx.norm <- cpm(atac.mx.raw,lib.size = atac.libsize)
      atac.dge <- DGEList(atac.mx.norm,group=group.names,lib.size = atac.libsize)
    }
    colnames(atac.mx.norm) <- colnames(atac.mx.raw)
    
    if (use.sva==T) {
      cat("Computing surrogate variables for data correction...\n")
      if (n.sva=="auto") {
        n.sv <- num.sv(cpm(atac.mx.raw,lib.size = atac.libsize)*1e6,atac.design)
      } else {
        n.sv <- as.integer(n.sva)
      }
      if (n.sv>0) {
        sv <- svaseq(cpm(atac.mx.raw,lib.size = atac.libsize,log = F)*1e6,atac.design,atac.design.null,n.sv = n.sv)$sv
        atac.design <- cbind(atac.design,sv)
        atac.design.null <- cbind(atac.design.null,sv)
        sv.data <- data.frame(sv,sampid=samp.names,libsize=atac.libsize,row.names = samp.names)
      } else {
        sv.data <- vector(mode = "numeric",length = 0L)
        cat("No significant surrogate variables found.\n")
      }
    } else {
      sv.data <- vector(mode = "numeric",length = 0L)
    }
    
    cat("\nEstimating dispersion parameters...\n")
    atac.dge <- estimateDisp(atac.dge, atac.design, robust=use.robust,tagwise = use.tagwise,trend.method = trend.method)
    
  } else {
    cat("Using pre-loaded DGE object, pass-thru. Omitting data load...\n")
    atac.mx.raw <- atac.dge$counts
    atac.design <- atac.dge$design
    
    if (norm.method!="cpm") {
      cat(paste("Applying",norm.method,"normalization...\n"))
      atac.mx.norm <- atac.mx.raw %*% diag(atac.dge$samples$norm.factors)
    } else {
      cat(paste("Applying",norm.method,"normalization...\n"))
      atac.mx.norm <- cpm(atac.mx.raw,lib.size = atac.libsize)
    }
    colnames(atac.mx.norm) <- colnames(atac.mx.raw)
    sv.data = NULL
  }
  
  cat("Fitting model...\n")
  atac.glm <- glmFit(atac.dge, atac.design)
  atac.glm.null <- glmFit(atac.dge, atac.design.null)
  
  if (!is.null(focus.contrast)) {
    cat("Executing GLM...\n")
    if (is.character(focus.contrast) | length(focus.contrast)==1) {
      atac.glmfit <- glmLRT(atac.glm, coef=focus.contrast)
    } else {
      atac.glmfit <- glmLRT(atac.glm, contrast=focus.contrast)
    }
    atac.glmtop <- topTags(atac.glmfit,n=nrow(atac.dge$counts),sort.by = "none", adjust.method = "fdr")$table
  } else {
    cat("No GLM test requested.\n")
    atac.glmtop = NULL
  }
  
  cat("Adjusting data...\n")
  prct = 0
  if (min(atac.mx.raw)<=0) prct = 1
  raw.data <- cpm(atac.mx.raw,prior.count = prct, log = F,lib.size = atac.libsize)
  
  norm.data <- log2(prct+atac.mx.norm)
  norm_mean <- matrix(rep(rowMeans(norm.data),ncol(norm.data)),nrow = nrow(norm.data))
  pred.data <- log2(prct+atac.glm$fitted.values) - log2(prct+atac.glm.null$fitted.values) + norm_mean
  adj.data <- norm.data - log2(prct+atac.glm.null$fitted.values) + norm_mean
  
  cat("Preparing output...\n")
  output <- list(raw=atac.mx.raw,
                 normd=norm.data,
                 pred=pred.data,
                 adj=adj.data,
                 glm=atac.glmtop,
                 name=run.name,
                 design=atac.design,
                 sva=sv.data,
                 dge=atac.dge,
                 glmfit=list(alt=atac.glm,null=atac.glm.null),
                 peaks=paste(atac.data.raw$chr,atac.data.raw$start,atac.data.raw$end,sep = "_"))
  cat("Done.\n")
  return(output)
}