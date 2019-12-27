browser_shot <- function(gene_query=NULL,
                         query_path=NULL,
                         bedgraph_dir=NULL,
                         ref_flank=2.5e5,
                         nref_peaks=1,
                         refpeak_fdr.cutoff=0.05,
                         atac.anno=NULL,
                         glm.results=list(Females=NULL,Males=NULL),
                         glm.dge=NULL,
                         atac.adj=NULL,
                         rna.adj=NULL,
                         valid.samps=NULL,
                         samp.stats=NULL,
                         draw.plot=T,
                         gpar=list(scanby="tss",
                                   rangex=c(-10000,10000),
                                   group.order = c("HY","HM","HO"),
                                   plot.style = c("line","area"),
                                   plot.se=T,
                                   area.alpha = 0.5,
                                   thinning.bw = 10,
                                   chrtrack.nudge=0.075,
                                   arrow.nudge = 0,
                                   arrow.freq = 10,
                                   bpaxistitlehjust = -1,
                                   grobwunit = 3.5,
                                   grobhunit = 0.5,
                                   max.fp = 10),
                         query_peaks_focus=T,
                         include_allcs=F,
                         include_footprints=T,
                         output_dir=NULL,
                         bedtools_path="~/bedtools/bin") {
  
  query = read.table(paste0(query_path,".qry"),sep = "\t",header = F,quote = "",col.names = c("chr.ref","start.ref","end.ref","GeneName","Specificity"),stringsAsFactors = F) %>%
    distinct()
  if (!nrow(query)) {
    print("Query is empty. Stopping now.")
    return(NULL)
  }
  
  require(ggbio)
  require(Homo.sapiens)
  require(ggplot2)
  require(reshape2)
  require(tibble)
  require(dplyr)
  
  data(genesymbol, package = "biovizBase")
  
  query_id <- basename(query_path)
  target_dir = paste0(query_path,"/",gene_query,"/")
  query_fname <- paste0(query_path,"/",gene_query,"/",gene_query,"_adjusted.data_atac_BySex.RData")
  
  this_query <- query %>% 
    filter(GeneName==gene_query)
  ct = unique(this_query$Specificity)
  if (!include_allcs) {
    cs <- switch(ct,
                 "CD4T"=c("CD4_naive","CD4_memory","PBMC"),
                 "CD8T"=c("CD8_naive","CD8_memory","PBMC"),
                 "naive_T"=c("CD4_naive","CD8_naive","PBMC"),
                 "memory_T"=c("CD4_memory","CD8_memory","PBMC"),
                 "CD3"=c("CD4_naive","CD8_naive","CD8_naive","CD8_memory","PBMC"),
                 "Lymphoid"=c("CD4_naive","CD4_memory","CD8_naive","CD8_memory","CD19","CD56","PBMC"),
                 "Myeloid"=c("CD14","PBMC"),
                 "Adaptive"=c("CD4_naive","CD8_naive","CD4_memory","CD8_memory","CD19","PBMC"),
                 "Innate"=c("CD14","CD56","PBMC"),
                 "CD4_naive"=c(ct,"PBMC"),
                 "CD4_memory"=c(ct,"PBMC"),
                 "CD8_naive"=c(ct,"PBMC"),
                 "CD8_memory"=c(ct,"PBMC"),
                 "CD19"=c(ct,"PBMC"),
                 "CD56"=c(ct,"PBMC"),
                 "CD14"=c(ct,"PBMC"),
                 ct)
  } else {
    cs <- c("CD14","CD56","CD4_naive","CD8_naive","CD4_memory","CD8_memory","CD19","PBMC")
  }
  if (query_peaks_focus) {
    # these peaks are used to draw links to Tss (they're highly correlated)
    query_highlight = (this_query %>%
                         mutate(peakco=paste(chr.ref,start.ref,end.ref,sep="_")))$peakco
  } else {
    query_highlight <- NULL
  }
  cat("Now processing gene:",gene_query,"from",query_id,"with cell specificity:",ct,"\n")
  
  samp.stats <- droplevels(samp.stats[valid.samps,])
  atac.glm.design=data.frame(glm.dge$design,row.names = rownames(glm.dge$samples))[valid.samps,]
  atac.glm.normfactors=glm.dge$samples[valid.samps,]$norm.factors
  atac.adj <- atac.adj[,valid.samps]
  rna.adj <- rna.adj[,intersect(colnames(rna.adj),valid.samps)]
  
  if (!basename(query_fname) %in% dir(path = dirname(query_fname))) {
    cat("Adjusted data not available. Proceeding to calculate it...")
    require(sva)
    require(edgeR)
    
    # Gathers metadata for queried peaks
    query_refs <- this_query %>%
      mutate(peakco=paste(chr.ref,start.ref,end.ref,sep = "_"),
             start.ref=start.ref-ref_flank,
             end.ref=end.ref+ref_flank,
             peakref=paste(GeneName,chr.ref,start.ref,end.ref,sep = "_")) %>%
      distinct(GeneName,peakco,peakref) %>%
      inner_join(atac.anno %>% dplyr::select(-GeneName),by="peakco") %>%
      inner_join(glm.results$Females %>% rename_all(list(~sub("$","_females_hohy",.))) %>% rownames_to_column("peakco"),by="peakco") %>%
      inner_join(glm.results$Males %>% rename_all(list(~sub("$","_males_hohy",.))) %>% rownames_to_column("peakco"),by="peakco")
    
    # Builds a list of queried bedGraph files
    cat(".")
    query_ls <- dir(bedgraph_dir,pattern = "*.query.txt",full.names = T)
    query_valid <- query_ls[sub("VHY402","VHY402D0",sub("\\_.*","",basename(query_ls))) %in% valid.samps]
    query_samps <- sub("VHY402","VHY402D0",sub("\\_.*","",basename(query_valid)))
    # Reads query bedGraphs (if not already read)
    query_bg <- lapply(setNames(as.list(query_samps),query_samps),function(n) {
      n=sub("VHY402D0","VHY402_D0",n)
      cat(n,"\n")
      read.table(query_valid[grepl(n,query_valid)],header = F,sep = "\t",quote = "",col.names = c("chr","start","end","score","chr.ref","start.ref","end.ref","GeneName","Specificity")) %>%
        mutate(peakref=paste(GeneName,chr.ref,start.ref,end.ref,sep = "_")) %>%
        dplyr::select(chr,start,end,score,peakref)
    })
    target_bg = paste(gene_query,"pooled.bedGraphUnion",sep = "_")
    
    system(paste("cd",getwd(),"; [[ ! -d",dirname(target_dir),"]] && mkdir",dirname(target_dir)))
    system(paste("cd",getwd(),"; [[ ! -d",target_dir,"]] && mkdir",target_dir))
    
    this_query_refs <- query_refs %>%
      filter(GeneName==gene_query) %>%
      dplyr::select(c(grep("^chromHMM",colnames(.),invert = T),grep(paste(paste("chromHMMstate",cs,sep = "_"),collapse="|"),colnames(.)))) %>%
      dplyr::select(c(grep("specificity",colnames(.),invert = T),grep("_specificity_fine",colnames(.))))
    
    # Retrieves bedGraph scores annotated to target gene peaks
    cat(".")
    top_peaks=nref_peaks # number of peaks to select. Peaks are sorted by combined significance (PValue in males and females), and the top most significants are selected for plotting
    this_target_ref = this_query_refs %>% 
      arrange(PValue_females_hohy*PValue_males_hohy) %>%
      dplyr::slice(top_peaks)
    # In addition to retrieving the data, subset bedGraph files are created and stored in a subfolder named after the target gene
    target_bg.lst <- do.call(rbind,lapply(setNames(as.list(names(query_bg)),names(query_bg)),function(n) {
      cat(n,"\n")
      X <- query_bg[[n]]
      bg = subset(X,peakref==this_target_ref$peakref)
      write.table(bg[,c("chr","start","end","score")],file = paste0(target_dir,n,"_",gene_query,".bedGraph"),quote = F,row.names = F,col.names = F,sep = "\t")
      return(bg)
    })) %>%
      rownames_to_column("sampid") %>%
      mutate(sampid=sub("\\..*","",sampid))
    
    # Uses bedtools to create a pooled bedGraph matrix from the union of the target-gene bedGraphs
    cat(".")
    system(command = paste0(bedtools_path,"/unionBedGraphs -header -names $(find ",
                            target_dir,
                            " -maxdepth 1 -name '*.bedGraph' | sed 's/.*\\///g' | sed 's/\\_.*//g') -i $(find ",
                            target_dir,
                            " -maxdepth 1 -name '*.bedGraph') > ",
                            target_dir,target_bg),
           intern=F)
    
    # Reads in and melts the pooled matrix bedGraph data
    cat(".")
    target_bg_pooled <- read.table(paste0(target_dir,gene_query,"_pooled.bedGraphUnion"),header = T,sep = "\t",quote = "") %>%
      dplyr::select("chrom","start","end",valid.samps) %>%
      rename_all(list(~sub("chrom","chr",.))) %>%
      melt(.,id.vars=c("chr","start","end"),variable.name="sampid",value.name="score") %>%
      mutate(sampid=as.character(sampid))
    
    # Batch correction
    cat(".")
    samp.stats$loglibsize = log1p(samp.stats$libsize)
    atac.design <- with(samp.stats,model.matrix(~sex + loglibsize + sex:age.group))
    atac.design.null <- with(samp.stats,model.matrix(~loglibsize))
    target_bg_pooled.df <- dcast(target_bg_pooled,formula = chr+start+end~sampid,value.var = "score") %>%
      filter(rowSums(.[,valid.samps])>0)
    target_bg_pooled.mx <- as.matrix(target_bg_pooled.df[,valid.samps]) %*% diag(samp.stats$libsize/1e7)
    target_bg_combadj.mx <- expm1(ComBat(dat = log1p(target_bg_pooled.mx),batch=samp.stats$batchdate,mod = atac.design))
    target_bg_combadj.mx[target_bg_combadj.mx < 0] <- 0 # zero negative values
    target_bg_combadj.mx <- setNames(data.frame(target_bg_combadj.mx %*% diag(colSums(target_bg_pooled.mx)/colSums(target_bg_combadj.mx))),valid.samps) # scale to original library sizes
    
    # Applies SV values calculated from base analysis, updates model
    cat(".")
    sv <- atac.glm.design[,last(colnames(atac.glm.design))]
    atac.design <- cbind(atac.design,sv)
    atac.design.null <- cbind(atac.design.null,sv)
    
    target_bg_dge <- with(samp.stats,DGEList(target_bg_combadj.mx,group=age.group,lib.size = libsize))
    target_bg_dge$samples$norm.factors <- atac.glm.normfactors
    target_bg_dge <- estimateDisp(target_bg_dge, atac.design, robust=T,tagwise = F,trend.method = "locfit")
    
    target_bg_glm <- glmFit(target_bg_dge, atac.design)
    target_bg_glm.null <- glmFit(target_bg_dge, atac.design.null)
    
    target_bg_norm <- log1p(target_bg_dge$counts %*% diag(target_bg_dge$samples$norm.factors))
    colnames(target_bg_norm) <- colnames(target_bg_dge$counts)
    target_bg_norm_mean <- matrix(rep(rowMeans(target_bg_norm),ncol(target_bg_norm)),nrow = nrow(target_bg_norm))
    target_bg_adj <- target_bg_norm - log1p(target_bg_glm.null$fitted.values) + target_bg_norm_mean
    
    target_bg_adjusted <- melt(cbind(target_bg_pooled.df[,c("chr","start","end")],target_bg_adj),id.vars=c("chr","start","end"),variable.name="sampid",value.name="score") %>%
      mutate(sampid=as.character(sampid))
    
    # Age by sex averaging
    cat(".")
    target_bg_pooled_agebysex <- target_bg_adjusted %>%
      inner_join(samp.stats,by="sampid") %>%
      mutate(sex=ifelse(sex=="F","Females","Males")) %>%
      group_by(chr,start,end,age.group,sex) %>%
      summarize(mean.score=mean(score),median.score=median(score),sd.score=sd(score,na.rm = T),se.score=sd(score,na.rm = T)/sqrt(n()),max.score=max(score)) %>%
      ungroup() %>%
      mutate(midpoint=0.5*(end+start))
    save(list = c("target_bg_pooled_agebysex","target_bg_adjusted","gene_query","query_id"),file = paste0(target_dir,gene_query,"_adjusted.data_atac_BySex.RData"))
    cat("done\n")
  } else {
    cat("Adjusted data already calculated, loading...")
    target_bg_pooled <- read.table(paste0(target_dir,gene_query,"_pooled.bedGraphUnion"),header = T,sep = "\t",quote = "") %>%
      dplyr::select("chrom","start","end",valid.samps) %>%
      rename_all(list(~sub("chrom","chr",.))) %>%
      melt(.,id.vars=c("chr","start","end"),variable.name="sampid",value.name="score") %>%
      mutate(sampid=as.character(sampid))
    cat(".")
    pooled_agebysex <- load(paste0(target_dir,gene_query,"_adjusted.data_atac_BySex.RData"))
    cat("done\n")
  }
  
  if (draw.plot) {
    cat("Calculating chromatin accessibility plot")
    if (gpar$scanby=="tss") {
      startx = gpar$rangex[1]
      endx = gpar$rangex[2]
      cat(paste0("as windowed range [",ifelse(startx<0,startx,paste0("+",startx)),"bp,",ifelse(endx<0,endx,paste0("+",endx)),"bp]"),"relative to",gene_query,"Tss\n")
      target_boundary <- range(genesymbol[gene_query],ignore.strand=T)
      start(target_boundary) <- max(start(target_boundary)+startx,min(target_bg_pooled$start))
      end(target_boundary) <- min(end(target_boundary)+endx,max(target_bg_pooled$end))
    } else if (gpar$scanby=="window") {
      target_boundary <- GRanges(with(target_bg_pooled_agebysex,paste0(as.character(chr)[1],":",min(start),"-",max(end))))
      cat("of entire covered range:",paste0(as.character(target_boundary@seqnames),":",as.character(target_boundary@ranges)),"\n")
    }
    
    target_bg_pooled_agebysex_thinned <- target_bg_pooled_agebysex %>%
      filter(midpoint %in% sort(unique(target_bg_pooled_agebysex$midpoint))[seq(1,nrow(target_bg_pooled_agebysex),gpar$thinning.bw)]) %>%
      filter(age.group!="HM") %>%
      mutate(age.group=factor(age.group,levels = gpar$group.order))
    
    bg <- ggplot(target_bg_pooled_agebysex_thinned %>%
                   filter(start>=target_boundary@ranges@start,
                          end<=with(target_boundary@ranges,start+width-1)),
                 aes(midpoint,mean.score,color=age.group,fill=age.group,group=age.group))
    if (gpar$plot.style[1]=="line") {
      bg <- bg + geom_line(size=0.25,alpha=0.85,na.rm = T)
    } else if (gpar$plot.style[1]=="area") {
      bg <- bg + geom_area(position="identity",alpha=gpar$area.alpha)
    }
    if (gpar$plot.se) bg <- bg + geom_ribbon(aes(ymin=mean.score-se.score,ymax=mean.score+se.score),alpha=0.1,size=0.05)
    bg <- bg +
      scale_x_continuous(labels = comma,expand = c(0,0),limits = with(target_boundary@ranges,c(start,end))) +
      scale_y_continuous(breaks = function(lim) c(0,floor(max(lim)))) +
      scale_color_manual(values = c(colors_age,colors_hmm18,colors_sex),guide=F) +
      scale_fill_manual(values = colors_age,guide=F) +
      facet_wrap(~sex,ncol=1) +
      labs(x=as.character(target_boundary@seqnames@values),y="Mean standardized accessibility") +
      theme_bw(base_size = 10) +
      theme(strip.text = element_text(size=16,angle=0,hjust=0),
            strip.background = element_blank(),
            axis.line.x = element_line(size=0.25),
            axis.title.x = element_text(hjust=0),
            axis.ticks.x = element_line(),
            panel.grid = element_blank(),
            plot.margin = margin(b=0))
    
    # Gene model track
    wr <- GRanges(paste(as.character(target_boundary@seqnames@values),paste(ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range,collapse = "-"),sep = ":"))
    wh <- selectByRanges(Homo.sapiens,wr,"SYMBOL")
    wh$SYMBOL[wh$SYMBOL=="MCUB"]<-"CCDC109B" ## synon
    # wh <- range(genesymbol[gene_query],ignore.strand=F)
    wg <- sapply(setNames(wh$SYMBOL,as.character(wh$SYMBOL)), function(g) range(genesymbol[intersect(genesymbol$symbol,g)]))
    wg <- wg[sapply(wg,length)>0]
    
    gm <- ggplot() + 
      geom_alignment(Homo.sapiens,which=wh,gap.geom="segment",stat="reduce",fill="black",label.color=NA,na.rm=T)
    gml <- gm +
      geom_text(data=gm$layers[[8]]$data,aes(midpoint,0,label=.labels),vjust=0.1,fontface="italic",size=4,na.rm = T) +
      scale_x_continuous(labels = comma,expand = c(0,0),limits = ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range) +
      theme_bw() +
      theme(aspect.ratio = 0.1,
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color=NA,fill=NA))
    
    gml <- gml +
      geom_line(aes(x=ifelse(rep(wg[[gene_query]]@strand=="+",2),
                             c(wg[[gene_query]]@ranges@start,min(ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range[2],with(wg[[gene_query]]@ranges,start+width-1))),
                             c(max(ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range[1],wg[[gene_query]]@ranges@start),with(wg[[gene_query]]@ranges,start+width-1))),
                    y=1)) +
      geom_point(aes(x=seq(from=ifelse(wg[[gene_query]]@strand=="+",
                                       wg[[gene_query]]@ranges@start,
                                       max(ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range[1],wg[[gene_query]]@ranges@start)),
                           to=ifelse(wg[[gene_query]]@strand=="+",
                                     min(ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range[2],with(wg[[gene_query]]@ranges,start+width-1)),
                                     with(wg[[gene_query]]@ranges,start+width-1)),
                           length.out = gpar$arrow.freq)[-1],
                     y=1),
                 shape=ifelse(wg[[gene_query]]@strand=="+",">","<"),size=3,position=position_nudge(y=gpar$arrow.nudge))
    
    # Peaks in interval
    wp <- melt(atac.anno[findOverlaps(wr,GRanges(atac.anno %>% distinct(chr,start,end)))@to,
                          c("chr","start","end","peakco",paste("chromHMMstate",cs,sep = "_"))],
                measure.vars=c("start","end"),variable.name="pos",value.name="coord") %>%
      mutate(sex="Females") %>%
      inner_join(glm.results$Females %>% rownames_to_column("peakco"),
                 by="peakco") %>%
      rbind(data.frame(.) %>% 
              mutate(sex="Males") %>%
              dplyr::select(grep(paste(colnames(glm.results$Females),collapse = "|"),colnames(.),invert = T)) %>%
              inner_join(glm.results$Males %>% rownames_to_column("peakco"),
                         by="peakco")) %>%
      mutate(isDifferential=FDR<=refpeak_fdr.cutoff) %>%
      left_join(data.frame(.) %>% 
                  group_by(peakco,sex) %>% 
                  summarize(logfc=mean(abs(logFC))) %>% 
                  group_by(sex) %>%
                  mutate(fcrank=rank(-logfc)) %>% 
                  ungroup() %>% 
                  filter(fcrank==nref_peaks) %>%
                  mutate(isTop=T) %>%
                  dplyr::select(-logfc,-fcrank),
                by=c("peakco","sex")) %>%
      mutate(isTop=!is.na(isTop))

    # Draw peak location tracks, colored using ENCODE chromHMM annotations
    for (ch in cs) {
      i = which(cs==ch)
      j = i%%2+1
      bg <- bg + 
        geom_hline(data=wp,yintercept=-0.1-gpar$chrtrack.nudge*(i-1),color=c("gainsboro","dimgray")[j],alpha=0.25,size=3.5) +
        geom_line(data=wp,aes_string(x="coord",y=as.character(-0.1-gpar$chrtrack.nudge*(i-1)),group="peakco",color=paste("chromHMMstate",ch,sep = "_")),alpha=1,size=3.5,lineend="butt",inherit.aes = F,na.rm = T) +
        annotate(geom = "text",x = ggplot_build(bg)$layout$panel_scales_x[[1]]$range$range[1],y = -0.1-gpar$chrtrack.nudge*(i-1), label=ch, hjust=0,vjust=0.5,size=2,fontface=ifelse(i==length(cs),"italic","plain"))
    }
    
    # Adds a significance indicator: stars based on Pvalue after filtering for FDR threshold
    bg <- bg +
      geom_text(data=wp %>% 
                  group_by(peakco,sex) %>% 
                  summarize(FDR=unique(FDR),PValue=unique(PValue),midpoint=mean(coord)) %>%
                  ungroup() %>%
                  mutate(siglab=ifelse(FDR<=refpeak_fdr.cutoff,ifelse(PValue<0.001,"***",ifelse(PValue<0.01,"**",ifelse(PValue<0.05,"*",""))),"")),
                aes(x=midpoint,y=0,label=siglab),inherit.aes = F,hjust=0.5,vjust=0.5)
    
    # Add a curved line connecting focus peak (eg peaks with accessibility highly correlated to expression), if using, to TSS (mind strand!)
    # if a peak is "differential" (as per FDR evaluation), solid curves are drawn to connect it to Tss
    # if a correlated peaks is not differential, dashed curves are drawn instead
    # Both diff and non-diff focus peaks are used in box plots. If no focus peaks are used, all peaks in the interval are used for box plots
    tss_loc <- ifelse(wg[[gene_query]]@strand=="+",wg[[gene_query]]@ranges@start,wg[[gene_query]]@ranges@start+wg[[gene_query]]@ranges@width-1)
    if (!is.null(query_highlight)) {
      wp2tss <- wp %>%
        filter(peakco %in% query_highlight) %>%
        group_by(sex,peakco) %>%
        mutate(coord=ifelse(pos=="start",mean(coord),tss_loc),
               isDifferential=all(isDifferential),
               isTop=all(isTop),
               logFC=mean(logFC),
               FDR=mean(FDR),
               PValue=mean(PValue)) %>%
        ungroup() %>%
        dcast(sex+peakco+isDifferential+isTop+logFC+FDR+PValue~pos,value.var = "coord") %>%
        mutate(
          lt=factor(ifelse(isDifferential,1,2)),
          fcscaled=scales::rescale(abs(logFC),to = c(0.1,1)))
      bg <- bg + 
        geom_curve(data=wp2tss,
                   aes(x=start,xend=end,y=-0.01,yend=-0.01,linetype=lt,size=fcscaled,color=sex),
                   # color="purple",
                   alpha=0.65,
                   lineend="round",
                   curvature=0.1,
                   angle=90,
                   ncp=3,
                   inherit.aes = F,
                   na.rm = T) +
        scale_linetype_manual(values = setNames(as.numeric(levels(wp2tss$lt)),levels(wp2tss$lt)),guide = F) +
        scale_size_identity(guide = F)
    }
    
    # Footprints in interval peaks, if using, by sex and age group (footprints must be precalculated for the query using findfootprints_at_queries.sh)
    # Only for TSS-style shots
    if (include_footprints & gpar$scanby=="tss") {
      fpfiles <- dir(path = paste0(query_path,"/motifs"),pattern = "*footprints.bed",full.names = T)
      fpsamps <- intersect(sub("\\_.*","",basename(fpfiles)),valid.samps)
      fptable <- do.call(rbind,lapply(fpsamps, function(n) {
        print(n)
        fpdata <- read.table(fpfiles[grep(n,fpfiles)],header = F,sep = "\t",quote = "",col.names = c("chr_fp","start_fp","end_fp","MotifName","Purity","strand_fp","chr","start","end","Target","ct","overlap")) %>%
          mutate(peakco=paste(chr,start,end,sep = "_")) %>%
          # filter(peakco %in% query_highlight,Target==gene_query) %>%
          filter(Target==gene_query) %>%
          mutate(MotifName=sub("\\.RC$","",sub("ver.*","",sub("Homer","",sub(".*full","",sub(".*dbd","",sub("MA0....","",MotifName)))))),
                 sampid=sub("VHY402","VHY402D0",n),
                 coord=0.5*(start+end))
      })) %>%
        inner_join(samp.stats,by="sampid") %>%
        filter(!is.na(age.group)) %>%
        distinct(Target,peakco,coord,MotifName,age.group,sex,sampid) %>%
        group_by(Target,coord,MotifName,age.group,sex) %>%
        summarize(n=n()) %>%
        ungroup() %>%
        dcast(Target+coord+MotifName+age.group~sex,value.var = "n",fill = 0) %>%
        melt(.,measure.vars=c("F","M"),variable.name="sex",value.name="n") %>%
        droplevels() %>%
        mutate(age.group=factor(age.group,levels=c("HY","HO"))) %>%
        dcast(Target+coord+MotifName+sex~age.group,value.var = "n",fill = 0,drop = F) %>%
        melt(.,measure.vars=c("HO","HY"),variable.name="age.group",value.name="n") %>%
        mutate(sex=ifelse(sex=="F","Females","Males"),
               MotifName=factor(MotifName),
               age.group=relevel(factor(age.group),ref = "HY")) %>%
        # this filters TFs to include only those with "sufficient" differnce between HO and HY. This will be done with a test (BiFET?) but for now, a cutoff where the top max.fp FPs by age difference are included. For completeness sake, if age-related differences of additional FPs above the top max.fp are equal (=same difference), they are all included 
        inner_join(data.frame(.) %>% 
                     group_by(MotifName,sex,age.group) %>%
                     mutate(p=mean(as.logical(n))) %>%
                     group_by(coord,MotifName,sex) %>%
                     summarize(d=abs(diff(n)),
                               dp=abs(diff(p))) %>% 
                     group_by(MotifName) %>%
                     summarize(D=max(d),Dp=max(dp)) %>%
                     ungroup() %>%
                     filter(D>0) %>%
                     arrange(-D,-Dp) %>%
                     filter(D>=(.[gpar$max.fp,])$D),
                   by="MotifName") %>%
        droplevels()
      
      bg <- bg +
        geom_point(data=fptable,aes(coord,-(as.numeric(MotifName)/10)-0.2*(length(cs)+1),alpha=n,shape=age.group),fill="black",color="black",position=position_dodge(width=0.025*target_boundary@ranges@width),stroke=0.25,na.rm = T) +
        geom_text(data=fptable,aes(coord+0.02*target_boundary@ranges@width,-(as.numeric(MotifName)/10)-0.2*(length(cs)+1),label=MotifName),hjust=0,size=1.5,color="black",fontface="italic",na.rm = T) +
        scale_shape_manual(values = c(25,21),guide=guide_legend(title = "Age group",title.theme = element_text(size=7),label.theme = element_text(size=4),override.aes = list(size=1))) +
        scale_alpha_continuous(range = c(0.1,0.9),guide=guide_legend(title = "N footprints",title.theme = element_text(size=7),label.theme = element_text(size=4),override.aes = list(size=1))) +
        geom_segment(data=fptable %>% group_by(coord) %>% summarize(mean_coord=mean(coord)) %>% ungroup(),
                     aes(x=mean_coord,
                         xend=mean_coord),
                     y=-0.21*(length(cs)+1),
                     yend=-length(levels(fptable$MotifName)),
                     alpha=0.25, color="dimgray", size=0.15,inherit.aes=F,na.rm=T) +
        theme(legend.position = c(0.01,0.01),
              legend.justification = c("left","bottom"))
    }
    
    # Peak accessibility box plots by sex
    if (!is.null(query_highlight)) {
      peaksel <- wp %>% 
        filter(peakco %in% wp2tss$peakco) %>%
        dplyr::select(sex,peakco) %>%
        distinct()
    } else {
      peaksel <- wp %>%
        dplyr::select(sex,peakco) %>%
        distinct()
    }
    target_bg_accessibility_bysex <- atac.adj %>%
      rownames_to_column("peakco") %>%
      melt(.,id.vars="peakco",variable.name="sampid",value.name="adj.atac") %>%
      mutate(sampid=as.character(sampid)) %>%
      inner_join(samp.stats,by="sampid") %>%
      mutate(sex=ifelse(sex=="F","Females","Males")) %>%
      inner_join(peaksel,by=c("sex","peakco"))

    abp <- ggplot(target_bg_accessibility_bysex,
                  aes(age.group,adj.atac,fill=age.group,color=age.group)) + 
      geom_boxplot(size=0.25,outlier.size = 0.5) + 
      facet_wrap(~sex,ncol=1) +
      scale_fill_manual(values = colors_age,guide=F) +
      scale_color_manual(values = colors_age_darken2,guide=F) +
      scale_y_continuous(position = "right") +
      labs(x="Age group",y="Adjusted peak accessibility") +
      theme_minimal(base_size = 10) +
      theme(axis.title.x = element_text(hjust=gpar$bpaxistitlehjust),
            axis.text.x = element_text(size=6,angle = 0),
            axis.line.x = element_line(size=0.25),
            plot.background = element_rect(fill="white",color="white"),
            strip.background = element_blank(),
            strip.text = element_text(size=16,angle = 0,color="white"),
            aspect.ratio = 3)
    
    # Gene expression box plots by sex
    target_bg_expression_bysex <- rna.adj[genesymbol[gene_query]$ensembl_id,] %>%
      rownames_to_column("EnsemblID") %>%
      melt(.,id.vars="EnsemblID",variable.name="sampid",value.name="adj.rna") %>%
      mutate(sampid=as.character(sampid)) %>%
      inner_join(samp.stats,by="sampid")
    
    rbp <- ggplot(target_bg_expression_bysex %>% mutate(sex=ifelse(sex=="F","Females","Males")),
                  aes(age.group,adj.rna,fill=age.group,color=age.group)) + 
      geom_boxplot(size=0.25,outlier.size = 0.5) + 
      facet_wrap(~sex,ncol=1) +
      scale_fill_manual(values = colors_age,guide=F) +
      scale_color_manual(values = colors_age_darken2,guide=F) +
      scale_y_continuous(position = "right") +
      labs(x=NULL,y=paste("Adjusted",gene_query,"expression")) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(size=6,angle = 0),
            axis.line.x = element_line(size=0.25),
            plot.background = element_rect(fill="white",color="white"),
            strip.background = element_blank(),
            strip.text = element_text(size=16,angle = 0,color="white"),
            aspect.ratio = 3)
    
    # Combo plot
    gtrack <- gtable_rbind(
      gtable_cbind(ggplotGrob(bg),
                   ggplotGrob(abp),
                   ggplotGrob(rbp)),
      gtable_cbind(ggplotGrob(gml),
                   ggplotGrob(ggplot(data=gm$layers[[8]]$data)+theme_bw()+theme(panel.border = element_rect(color=NA,fill=NA))),
                   ggplotGrob(ggplot(data=gm$layers[[8]]$data)+theme_bw()+theme(panel.border = element_rect(color=NA,fill=NA)))))
    gtrack$widths[c(14,23)] <- unit(gpar$grobwunit,units = "grobwidth",data = ggplotGrob(bg))
    gtrack$heights[25] <- unit(gpar$grobhunit,units = "grobheight",data = ggplotGrob(bg))
    
    # Export PDF
    if (is.null(output_dir)) output_dir <- paste0(query_path,"/",gene_query,"/")
    pdf(paste0(output_dir,gene_query,"_",gpar$scanby,".shot",ifelse(gpar$scanby=="window","",paste0(".",startx,".",endx)),".pdf"))
      plot(gtrack)
    dev.off()
    return(gtrack)
  }
}
