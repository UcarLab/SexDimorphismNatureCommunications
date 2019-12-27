# Function to run differential analysis (DA) on RNA-seq for the sex:age data using the specified models and parameters.
# Function assumes RNAseq data objects obtained after benchmarking are loaded in the working environment.
# Example call used in these analyses:
# 
# # Global parameters used in DA
# global.pars <- list(use.adjustment="libsize", # no|libsize|frips|depth
#                     use.combat=T,
#                     batch.variable="batchdate",
#                     use.sva=T,
#                     n.sva="auto",
#                     use.robust=T,
#                     use.tagwise=T,
#                     norm.method="TMM",
#                     trend.method="locfit",
#                     min.expressed=2,
#                     protein_coding_only=T,
#                     renormalize=T,
#                     run.name=run.name)
# if (!global.pars$use.sva) global.pars$n.sva=0
# 
# design <- paste("model.matrix(~sex +",global.pars$batch.variable,"+ sex:age.group)")
# design.null <- paste0("model.matrix(~",global.pars$batch.variable,")")
# 
# # Data frame of relevant factors, individually loaded from filtered_global.stats_narrowPeaks_filtered.RData
# metadata.factors <- data.frame(sex=as.character(sex),
#                                age=age,
#                                age.group=as.character(age.group),
#                                batchdate=as.character(batchdate),
#                                season=season,
#                                sampid=as.character(samps),
#                                libsize=log2(libsize.rna.pbmc),
#                                stringsAsFactors = F)
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
# rna.design <- with(metadata.factors,eval(parse(text=design)))
# rna.design.null <- with(metadata.factors,eval(parse(text=design.null)))
# 
# # Batch correction via ComBat
# transcripts <- rna.pbmc
# if (globalPars$use.combat) {
#   cat("Applying batch correction via ComBat, batch variable:",globalPars$batch.variable)
#   library(sva)
#   transcripts.combadj <- expm1(ComBat(dat = as.matrix(log1p(transcripts[,metadata.factors$sampid])),batch = get(globalPars$batch.variable)[metadata.factors$sampid],mod = rna.design))
#   transcripts.combadj[transcripts.combadj < 0] <- 0 # zero negative values
#   transcripts <-setNames(data.frame(as.matrix(transcripts.combadj) %*% diag(libsize/colSums(transcripts.combadj))),names(libsize)) # scale to original library sizes
# }
# 
# # Definition of contrasts to test, variable whether DA compare age groups by sex, or sexes by age group
# # This is for analysis by sex
# focus.contrast.factors = "age"
# focus.contrast.set <- c(as.list(colnames(rna.design)[grep(paste(focus.contrast.factors,collapse="|"),colnames(rna.design))]),
#                         list(c(as.numeric(grepl("sexF\\:",colnames(rna.design))),rep(0,ncol(da.list$design)-ncol(rna.design))),
#                              c(as.numeric(grepl("sexM\\:",colnames(rna.design))),rep(0,ncol(da.list$design)-ncol(rna.design)))))
# focus.contrast.set <- lapply(focus.contrast.set,function(x) {
#   if (is.numeric(x)) {
#     x[min(which(as.logical(x)))] <- -1
#   }
#   return(x)
# })
# focus.contrast.set.names <- sapply(focus.contrast.set, function(x) {
#   if (is.numeric(x)) {
#     n <- sub("sex.*age","Age",sub("^sexM\\:age.groupage","Males_Age",sub("^sexF\\:age.groupage","Females_Age",paste(rev(colnames(rna.design)[as.logical(x)]),collapse="x"))))
#   } else {
#     n <- paste0(sub("M","Males_Age",sub("F","Females_Age",sub("sex","",sub("\\_.*age","",sub("\\:","_",x))))),"xAge1")
#   }
# })
# 
# # Differential analysis routine call
# da_bysex <- lapply(setNames(focus.contrast.set,focus.contrast.set.names), function(focus.contrast) {
#   cat("Now processing:",focus.contrast,"...\n")
#   da.list <- da_expression(rna.data.raw=transcripts,
#                            rna.expressed=expressed,
#                            rna.anno=anno,
#                            rna.design=rna.design,
#                            rna.design.null=rna.design.null,
#                            focus.contrast=focus.contrast,
#                            samp.names=metadata.factors$sampid,
#                            run.name=globalPars$run.name,
#                            group.names=age.group,
#                            rna.libsize=libsize,
#                            use.adjustment=globalPars$use.adjustment,
#                            min.expressed=globalPars$min.expressed,
#                            protein_coding_only=globalPars$protein_coding_only,
#                            renormalize=globalPars$renormalize,
#                            use.sva=globalPars$use.sva,
#                            n.sva=globalPars$n.sva,
#                            use.robust=globalPars$use.robust,
#                            use.tagwise = globalPars$use.tagwise,
#                            norm.method = globalPars$norm.method,
#                            trend.method = globalPars$trend.method)
#   cat("Found",sum(da.list$glm$FDR<=fdr.thresh),"hits at",100*fdr.thresh,"% FDR\n")
#   da.output <- list(da.list=da.list,
#                     metadata.factors=metadata.factors,
#                     rna.design=rna.design,
#                     rna.design.null=rna.design.null,
#                     focus.contrast=focus.contrast,
#                     global.pars=global.pars,
#                     age.group=age.group)
#   return(da.output)
# })
#################################################################################################################################################################

da_expression <- function(rna.dge=NULL,
                          rna.data.raw=NULL,
                          rna.expressed=NULL,
                          rna.anno=NULL,
                          rna.design=NULL,
                          rna.design.null=NULL,
                          focus.contrast=NULL,
                          rna.libsize=NULL,
                          samp.names=NULL,
                          group.names=NULL,
                          run.name=NULL,
                          exclude.samples=NULL,
                          use.adjustment=c("no","libsize"),
                          min.expressed=1,
                          protein_coding_only=FALSE,
                          renormalize=TRUE,
                          use.sva=F,
                          n.sva="auto",
                          use.robust=T,
                          use.tagwise=F,
                          norm.method=c("TMM","RLE","upperquartile","cpm","none"),
                          trend.method="none") {
  
  require(edgeR)
  require(sva)
  require(reshape2)
  require(dplyr)
  
  protein_coding_cats = c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_M_gene", "IG_V_gene", "IG_Z_gene", "nonsense_mediated_decay", "nontranslating_CDS", "non_stop_decay", "polymorphic_pseudogene", "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", "TR_J_gene", "TR_V_gene")
  # (ref: http://useast.ensembl.org/Help/Glossary?id=275)
  
  if (is.null(rna.dge)) {
    if (!is.null(exclude.samples)) {
      samp.names <- setdiff(samp.names,exclude.samples)
    }
    rna.libsize <- rna.libsize[samp.names]
    group.names <- droplevels(group.names[samp.names])
    
    n_init <- nrow(rna.data.raw)
    rna.data.raw <- rna.data.raw[,samp.names] %>% 
      mutate(EnsemblID=rownames(.)) %>%
      select(EnsemblID,everything()) %>%
      left_join(rna.anno %>% select(EnsemblID,GeneType),by = "EnsemblID") %>%
      left_join(rna.expressed[,samp.names] %>% 
                  mutate(EnsemblID=rownames(.)) %>%
                  slice(match(EnsemblID,rownames(rna.data.raw))) %>%
                  mutate(Total.counts=apply(.[,samp.names],1,sum)) %>% 
                  select(EnsemblID,Total.counts),
                by = "EnsemblID") %>%
      filter(ifelse(rep(protein_coding_only,nrow(.)),GeneType %in% protein_coding_cats,TRUE),!Total.counts<min.expressed)
    if (renormalize) rna.libsize=colSums(rna.data.raw[,samp.names])
    
    rna.mx.raw <- as.matrix(rna.data.raw[,samp.names])
    rownames(rna.mx.raw) <- rna.data.raw$EnsemblID
    
    cat(paste("Run name:",run.name,"\n"))
    cat(paste("Input data:",deparse(quote(rna.data.raw)),"\n"))
    cat(paste("Sample names:",paste(samp.names,collapse = "-"),"\n"))
    cat(paste("group names:",paste(levels(group.names),collapse = ","),"\n"))
    cat(paste("Filter 1: each transcript must be expressed in at least",min.expressed,"sample(s) to enter the analysis.\n"))
    cat(paste("Filter 2: analysis based on",ifelse(protein_coding_only,"protein-coding","all"),"genes.\n"))
    cat(paste("==>",nrow(rna.data.raw),"out of",n_init,"transcripts entered the analysis.\n"))
    cat(paste("Adjusted: ",use.adjustment,"\n"))
    cat(paste("Normalization method: ",norm.method,"\n"))
    cat(paste("Using Tagwise dispersion estimates:",use.tagwise,"\n"))
    cat(paste("Dispersion estimation trend method:",trend.method,"\n"))
    cat(paste("Using SVA:",use.sva,"\n"))
    if (use.sva) cat(paste("N SV:",n.sva,"\n"))
    cat(paste("Using robust estimate of dispersion factors:",use.robust,"\n"))
    
    ## edgeR DA of joint cohort RNAseq genes/transcripts
    cat("Building data matrix...\n")
    if (norm.method!="cpm") {
      rna.dge <- DGEList(rna.mx.raw,group=group.names,lib.size = rna.libsize)
      cat(paste("Applying",norm.method,"normalization...\n"))
      rna.dge <- calcNormFactors(rna.dge,method = norm.method);
      rna.mx.norm <- rna.mx.raw %*% diag(rna.dge$samples$norm.factors)
    } else {
      cat(paste("Applying",norm.method,"normalization...\n"))
      rna.mx.norm <- cpm(rna.mx.raw,lib.size = rna.libsize)
      rna.dge <- DGEList(rna.mx.norm,group=group.names,lib.size = rna.libsize)
    }
    colnames(rna.mx.norm) <- colnames(rna.mx.raw)
    
    if (use.sva==T) {
      cat("Computing surrogate variables for data correction...\n")
      if (n.sva=="auto") {
        n.sv <- num.sv(cpm(rna.mx.raw,lib.size = rna.libsize)*1e6,rna.design)
      } else {
        n.sv <- as.integer(n.sva)
      }
      if (n.sv>0) {
        sv <- svaseq(cpm(rna.mx.raw,lib.size = rna.libsize,log = F)*1e6,rna.design,rna.design.null,n.sv = n.sv)$sv
        rna.design <- cbind(rna.design,sv)
        rna.design.null <- cbind(rna.design.null,sv)
        sv.data <- data.frame(sv,sampid=samp.names,libsize=rna.libsize,row.names = samp.names)
      } else {
        sv.data <- vector(mode = "numeric",length = 0L)
        cat("No significant surrogate variables found.\n")
      }
    } else {
      sv.data <- vector(mode = "numeric",length = 0L)
    }
    
    cat("\nEstimating dispersion parameters...\n")
    rna.dge <- estimateDisp(rna.dge, rna.design, robust=use.robust,tagwise = use.tagwise,trend.method = trend.method)
    
  } else {
    cat("Using pre-loaded DGE object, pass-thru. Omitting data load...\n")
    rna.mx.raw <- rna.dge$counts
    rna.design <- rna.dge$design
    
    if (norm.method!="cpm") {
      cat(paste("Applying",norm.method,"normalization...\n"))
      rna.mx.norm <- rna.mx.raw %*% diag(rna.dge$samples$norm.factors)
    } else {
      cat(paste("Applying",norm.method,"normalization...\n"))
      rna.mx.norm <- cpm(rna.mx.raw,lib.size = rna.libsize)
    }
    colnames(rna.mx.norm) <- colnames(rna.mx.raw)
    sv.data = NULL
  }
  
  cat("Fitting model...\n")
  rna.glm <- glmFit(rna.dge, rna.design)
  rna.glm.null <- glmFit(rna.dge, rna.design.null)
  
  # Limma code example
  # rna.glm <- lmFit(log2(rna.dge$counts), rna.design)
  # rna.confit <- contrasts.fit(rna.glm,coef=focus.contrast)
  # rna.glmfit <- eBayes(rna.glm)
  # rna.glmtop <- topTable(rna.glmfit,coef = focus.contrast,number = Inf,adjust.method = "fdr",sort.by = "none")
  
  if (!is.null(focus.contrast)) {
    cat("Executing GLM...\n")
    if (is.character(focus.contrast) | length(focus.contrast)==1) {
      rna.glmfit <- glmLRT(rna.glm, coef=focus.contrast)
    } else {
      rna.glmfit <- glmLRT(rna.glm, contrast=focus.contrast)
    }
    rna.glmtop <- topTags(rna.glmfit,n=nrow(rna.dge$counts),sort.by = "none", adjust.method = "fdr")$table
  } else {
    cat("No GLM test requested.\n")
    rna.glmtop = NULL
  }
  
  cat("Adjusting data...\n")
  prct = 0
  if (min(rna.mx.raw)<=0) prct = 1
  raw.data <- cpm(rna.mx.raw,prior.count = prct, log = F,lib.size = rna.libsize)
  
  norm.data <- log2(prct+rna.mx.norm)
  norm_mean <- matrix(rep(rowMeans(norm.data),ncol(norm.data)),nrow = nrow(norm.data))
  pred.data <- log2(prct+rna.glmfit$fitted.values) - log2(prct+rna.glm.null$fitted.values) + norm_mean
  adj.data <- norm.data - log2(prct+rna.glm.null$fitted.values) + norm_mean
  
  cat("Preparing output...\n")
  output <- list(raw=raw.data,
                 normd=norm.data,
                 pred=pred.data,
                 adj=adj.data,
                 glm=rna.glmtop,
                 name=run.name,
                 design=rna.design,
                 sva=sv.data,
                 dge=rna.dge,
                 glmfit=list(alt=rna.glm,null=rna.glm.null),
                 transcripts=rownames(rna.data.raw))
  cat("Done.\n")
  return(output)
}