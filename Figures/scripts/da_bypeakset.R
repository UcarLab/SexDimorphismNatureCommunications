# Function to run differential analysis on ATAC-seq for the sex:age data using the specified models and parameters

da_bypeakset <- function(atac.dge=NULL,
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
                         norm.method=c("TMM","RLE","upperquartile","cpm","none"),
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
      # left_join(atac.peakcalls %>% select(chr,start,end,Total.counts),by = c("chr","start","end")) %>%
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
    
    ## edgeR DA of joint cohort ATAC peaks
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
        # sv=fastICA(sv,n.comp = 3)$S
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
  
  # Limma code example
  # atac.glm <- lmFit(log2(atac.dge$counts), atac.design)
  # atac.confit <- contrasts.fit(atac.glm,coef=focus.contrast)
  # atac.glmfit <- eBayes(atac.glm)
  # atac.glmtop <- topTable(atac.glmfit,coef = focus.contrast,number = Inf,adjust.method = "fdr",sort.by = "none")
  
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