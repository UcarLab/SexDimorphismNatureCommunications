### export source tables
# Fig 1 (load Figure1.Rmd files)
##############################################################################################################################
Fig1b <- with(atacseq,
              Vall$pcstats %>% 
                filter(grepl("PC1",PCn))) %>%
  select(sampid,sex,age.group,PCn,score)
savesource("Fig1b")
##############################################################################################################################
Fig1c <- gs.combo %>%
  filter(geneset=="vp2008") %>%
  dcast(., Module.Name~data.source+direction,value.var="signed.hypergeom.p") %>%
  filter_all(function(x) !is.na(x)) %>%
  melt(.,id.vars=c("Module.Name"),variable.name = "data.source",value.name = "signed.hypergeom.p") %>%
  mutate(direction=sub(".*\\_","",data.source),
         data.source=sub("\\_.*","",data.source)) %>%
  arrange(signed.hypergeom.p) %>% 
  mutate(Module.Name=factor(Module.Name,levels = unique(.$Module.Name)),
         signed.hypergeom.p=ifelse(abs(signed.hypergeom.p)>40,sign(signed.hypergeom.p)*40,signed.hypergeom.p))
# Note: in fig1c and figs1c the scale limits were cut at loginvp=40 for visibility. 
# Pvalues larger than this were set at 40 (1c) and 20 (s1c). Therefore, the scale should be written as something like "40+".
# Also, no negative pvalues should be used. Instead, clarify that the red/blue colors refer to the direction of 
# enrichment, not the Pvalue signs.
# Also, notice the final figure excludes the modules that were not significant for either ATAC or RNAseq 
# (significance threshold=5%FDR)
savesource("Fig1c")
##############################################################################################################################
FigS1a <- rnaseq$Vall$pcdata %>% select(sampid,sex,age.group,everything(),-batchdate,-nreps)
savesource("FigS1a")
FigS1b <- atacseq$Vall$pcdata %>% select(sampid,sex,age.group,everything(),-batchdate,-nreps)
# Note: these tables include sample ids (non-anonimized and batch info that may not be relefvant for this paper)
savesource("FigS1b")
##############################################################################################################################
FigS1c <- gs.combo %>%
  filter(geneset=="scrnaseq_pbmc_simple_specific") %>%
  dcast(., Module.Name~data.source+direction,value.var="signed.hypergeom.p") %>%
  filter_all(function(x) !is.na(x)) %>%
  melt(.,id.vars=c("Module.Name"),variable.name = "data.source",value.name = "signed.hypergeom.p") %>%
  mutate(direction=sub(".*\\_","",data.source),
         data.source=sub("\\_.*","",data.source)) %>%
  arrange(data.source,signed.hypergeom.p) %>% 
  mutate(Module.Name=factor(Module.Name,levels = unique(.$Module.Name)),
         signed.hypergeom.p=ifelse(abs(signed.hypergeom.p)>20,sign(signed.hypergeom.p)*20,signed.hypergeom.p))
savesource("FigS1c")
##############################################################################################################################
FigS1g1 <- age.dist %>% select(-ageLabel)
savesource("FigS1g1")
FigS1g2 <- frailty %>% select(-flab)
savesource("FigS1g2")
##############################################################################################################################


# Fig 2 (load Figure2.Rmd and Figure3.Rmd files)
##############################################################################################################################
Fig2a <- da.glm_agebysex %>% 
  mutate(Contrast_label=factor(sub("_",", ",sub("x"," vs ",sub("Age1","HY",sub("Age2","HM",sub("Age3","HO",Contrast_label))))))) %>%
  select(-nhits,-Contrast_color,-Contrast_label)
savesource("Fig2a")
##############################################################################################################################
Fig2b <- da.chromhmm.dist_agebysex %>%
  select(peakco,everything,-label,-nhits,-label_brief,-n) %>%
  mutate(direction=sub(" peaks in","",direction))
savesource("Fig2b")
##############################################################################################################################
Fig2c_S2c <- da.chromhmm.spec_agebysex %>% 
  select(peakco,SexContrast,direction,paste0("RE_specificity_",ctset)) %>%
  rename_all(funs(sub("RE_specificity.*","Specificity",.))) %>%
  filter(!is.na(Specificity)) %>%
  mutate(Specificity=ifelse(grepl("^B|Plasma",Specificity),"CD19",Specificity),
         Specificity=factor(ifelse(Specificity %in% c(focus.celltypes),Specificity,"Other")),
         direction=sub(" peaks in","",direction))
savesource("Fig2c_S2c")
##############################################################################################################################
load("H:/Downloads/IL7R_adjusted.data_atac_BySex.RData")
Fig2d <- target_bg_adjusted %>%
  mutate(gene.focus=gene_query)
# Replace sample names for anonymized versions
savesource("Fig2d")
##############################################################################################################################
# @fig3
Fig2e <- combo.tested.modules.results[[1]] %>%
  filter(Module.Name %in% c(focus.celltypes,"Monocytes","NK_cells","Bcells","acCD8_Tcells","Tcells","Naive_Tcells","Cytotoxic cells","B cells","Inflammation I","Inflammation II","T cells","Myeloid lineage 1","Myeloid lineage 2","B_naive")) %>%
  filter(grepl("Age3xAge1",Contrast_label)) %>%
  mutate(Contrast_color=ifelse(signed.hypergeom.p<0,ifelse(grepl("Age1",.$Contrast_label),colors_age["HY"],colors_age["HM"]),ifelse(grepl("Age3",.$Contrast_label),colors_age["HO"],colors_age["HM"])),
         Contrast_alpha=ifelse(hypergeom.p>combo.threshP,1,0.5),
         Contrast_label=factor(sub("\nAge3xAge1","",Contrast_label)),
         Module.Label=sub(" genes","",as.character(Module.Label))) %>%
  group_by(Module.Label) %>%
  mutate(m=mean(signed.hypergeom.p)) %>%
  ungroup() %>%
  arrange(m) %>%
  mutate(Module.Label=factor(Module.Label,levels = unique(as.character(.$Module.Label)))) %>% 
  left_join(data.frame(.) %>%
              group_by(Module.Label) %>%
              summarize(m=mean(m)) %>%
              mutate(indx=row_number()) %>%
              select(Module.Label,indx),
            by="Module.Label") %>% 
  ungroup() %>% 
  filter(indx<=ntop.gs | indx>=max(indx)-ntop.gs) %>% 
  mutate(signed.hypergeom.p = ifelse(hypergeom.p>maxlogP,sign(signed.hypergeom.p)*maxlogP,signed.hypergeom.p)) %>%
  filter(data.source=="RNA") %>%
  select(celltype="Module.Name",direction="gene.set",gene.count,Total.Freq,hypergeom.p,hypergeom.fdr,signed.hypergeom.p)
savesource("Fig2e")
##############################################################################################################################
# @fig3
gs="T Cells"
transcription <- da.glm_agebysex.rna %>% filter(FDR<=fdr.thresh.rna) %>% inner_join(anno.rna.pbmc,by="EnsemblID")
commonGenes = intersect(gsets_gs[[gs]]$GeneName,transcription$GeneName)
topGSdata <- transcription %>% 
  filter(GeneName %in% commonGenes) %>%
  arrange(GeneName,PValue) %>%
  filter(!duplicated(GeneName)) %>%
  merge(adj.data.rna %>% select(EnsemblID,rnaseq$samps),by="EnsemblID")
Fig2f <- do.call(rbind,lapply(list(Females="F",Males="M"), function(sx) {
  dt <- topGSdata %>% 
    select(GeneName,rnaseq$samps[rnaseq$sex==sx]) %>% 
    column_to_rownames("GeneName") %>% 
    t() %>%
    data.frame() %>%
    rownames_to_column("sampid") %>%
    inner_join(data.frame(age.group=rnaseq$age.group,sex=rnaseq$sex) %>% rownames_to_column("sampid"),by="sampid") %>%
    select(sampid,sex,age.group,everything())
})) %>%
  data.frame(.,row.names = NULL)
# Note: replace sample names with anonymized names
savesource("Fig2f")
##############################################################################################################################
FigS2a <- list(ATAC=da.glm_venn_bypeak,RNA=da.glm_venn_bytss)
# R list, saved as compressed R file
save(FigS2a,file = "source_data_final/Fig_S2a.RData")
##############################################################################################################################
FigS2b <- da.chromhmm.spec_agebysex %>% 
  select(peakco,SexContrast,direction,paste0("RE_specificity_",ctset)) %>%
  rename_all(funs(sub("RE_specificity.*","Specificity",.))) %>%
  mutate(Specificity=factor(ifelse(is.na(Specificity),"Cell-agnostic state","Cell-specific state")),
         direction=sub(" peaks in","",direction))
savesource("FigS2b")
##############################################################################################################################
FigS2d <- da.glm_agebysex %>% 
  mutate(SexContrast=sub("\\_.*","",Contrast),
         direction=ifelse(logFC<0,"Closing peaks","Opening peaks")) %>%
  inner_join(atac.allpeaks.anno.pbmc %>% 
               select(peakco,paste0("RE_specificity_",ctset)) %>%
               rename_all(funs(sub("RE_specificity.*","Specificity",.))) %>%
               filter(!is.na(Specificity)),
             by="peakco") %>% 
  mutate(Specificity=ifelse(grepl("^B|Plasma",Specificity),"CD19",Specificity),
         Specificity=factor(ifelse(Specificity %in% c(focus.celltypes),Specificity,"Other")),
         direction=sub(" in$","",direction)) %>%
  group_by(SexContrast,Specificity) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  mutate(Specificity=factor(Specificity,levels = c(focus.celltypes,"Other")),
         label=factor(paste0(Specificity,"\nn = ",n))) %>%
  mutate(label=factor(label,levels = levels(label)[c(1,5,2,4,3,7,6,8)])) %>%
  filter(Specificity!="Other") %>%
  select(-Contrast,-Contrast_color,-nhits,-Contrast_label,-label,-n)
savesource("FigS2d")
##############################################################################################################################
FigS2e <- da.chromhmm_enrichment_agebysex %>%
  filter(celltype %in% focus.celltypes) %>%
  mutate(celltype=factor(celltype,levels = focus.celltypes),
         direction=paste(sub("with aging","in",direction),tolower(SexContrast)),
         hypergeom.p=ifelse(hypergeom.p > -log10(1e-50),-log10(1e-75),hypergeom.p),
         signed.hypergeom.p=sign(signed.hypergeom.p)*hypergeom.p) %>%
  select(-Contrast,-SexContrast,-celltype_set,-AgeContrast)
# Note: somewhere clarify that the signs represent which set of peaks (opening/closing) is being tested for enrichment
savesource("FigS2e")
  
# Fig3 (load Figure4.Rmd)
##############################################################################################################################
# @fig4
Fig3a_3b <- logfc_bysex %>%
  select(chr,start,end,GeneName,Annotation,DistancetoTSS,Specificity="Module.State",Females.atac,Males.atac,Females.rna,Males.rna)
savesource("Fig3a_3b")
##############################################################################################################################
# @fig4
ctset = "state"
assign(paste0("logfc_bysex_",ctset),logfc_bysex_byCS[[ctset]]$logfc_bysex %>%
         filter(Specificity.Match %in% cs.included[[ctset]]$Specificity.Alias) %>%
         mutate(Specificity=factor(Specificity,levels = c(cs.included[[ctset]]$Specificity,"Other")),
                Specificity.Match=factor(Specificity.Match,levels = cs.included[[ctset]]$Specificity.Alias)))
Fig3c <- logfc_bysex_state %>%
  filter(Specificity!="Other") %>%
  droplevels() %>%
  select(chr,start,end,GeneName,Annotation,DistancetoTSS,Specificity,Females.atac,Males.atac)
savesource("Fig3c")
##############################################################################################################################
# @fig4
gs = "scrnaseq_pbmc_simple_specific"
assign(paste0("logfc_bysex_",gs),logfc_bysex_byGS[[gs]]$logfc_bysex %>%
         left_join(gs.modules.included[[gs]],by="Module.Name") %>%
         filter(Module.Match %in% gs.modules.included[[gs]]$Module.Name) %>%
         mutate(Module.Name=factor(Module.Name,levels = c(gs.modules.included[[gs]]$Module.Name,"Other")),
                Module.Match=factor(Module.Match, levels = gs.modules.included[[gs]]$Module.Name)))
Fig3d <- logfc_bysex_scrnaseq_pbmc_simple_specific %>%
  filter(Module.Name!="Other") %>%
  droplevels() %>%
  select(GeneName,Module.Name,Females.rna,Males.rna)
savesource("Fig3d")
##############################################################################################################################
# @ fig4
Fig3e <- adj.data.joint.highlights %>%
  select(sampid,sex,age.group,DataSource,chr,start,end,GeneName,adj.score)
# Note: replace sample names with anonymized names
savesource("Fig3e")
##############################################################################################################################
load("H:/Downloads/S100A9_adjusted.data_atac_BySex.RData")
FigS3a <- target_bg_adjusted %>%
  mutate(gene.focus=gene_query)
# Replace sample names for anonymized versions
savesource("FigS3a")
##############################################################################################################################
load("H:/Downloads/CD79B_adjusted.data_atac_BySex.RData")
FigS3b <- target_bg_adjusted %>%
  mutate(gene.focus=gene_query)
# Replace sample names for anonymized versions
savesource("FigS3b")
##############################################################################################################################
load("H:/Downloads/GNLY_adjusted.data_atac_BySex.RData")
FigS3c <- target_bg_adjusted %>%
  mutate(gene.focus=gene_query)
# Replace sample names for anonymized versions
savesource("FigS3c")
##############################################################################################################################

# Fig 4 (load Figure3.Rmd and Figures5&6.Rmd files)
##############################################################################################################################
# @ fig5&6
Fig4af <- predicted_bysex$females
savesource("Fig4af")
Fig4am <- predicted_bysex$males
savesource("Fig4am")
##############################################################################################################################
# @ fig5&6
focus_cs <-c("Monocytes","CD14","NK_cells","CD56","CD4_naive","CD4_memory","CD8_naive","CD8_memory","CD3","B_naive","B_memory","Plasma_cells","CD19")
focus_cs_names <-c("Monocytes","Monocytes","NK cells","NK cells","CD4 naive","CD4 memory","CD8 naive","CD8 memory","T cells","B naive","B memory","Plasma cells","B cells")
Fig4b <- rbind(chrspecific.enrichment_bysex.atac$females$sublineage$all %>%
                 select(sex,everything(),-scope),
               chrspecific.enrichment_bysex.atac$males$sublineage$all %>%
                 select(sex,everything(),-scope)) %>%
  filter(celltype %in% focus_cs) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_cs,to = focus_cs_names))
savesource("Fig4b")
##############################################################################################################################
# @ fig5&6
Fig4cf <- predicted_bysex.rna$females
savesource("Fig4cf")
Fig4cm <- predicted_bysex.rna$males
savesource("Fig4cm")
##############################################################################################################################
# @ fig5&6
focus_ct <-c("Monocytes","NK_cells","Naive_Tcells","CD4_naive","CD8_naive","Tcells","acCD8_Tcells","Bcells","B_naive","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","CD4 naive","CD8 naive","T cells","Act CD8 T cells","B cells","B cells","Plasma cells")
Fig4d <- rbind(geneset.enrichment_bysex.rna$females$scrnaseq_pbmc_simple_specific$all %>%
                 select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set"),
               geneset.enrichment_bysex.rna$males$scrnaseq_pbmc_simple_specific$all %>%
                 select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set")) %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))
savesource("Fig4d")
##############################################################################################################################
# @ fig5&6 - breakpoints_atac chunk
Fig4e_S6c <- do.call(rbind,lapply(c("females","males"), function(sx) {
  X <- allBreaks[[sx]]
  do.call(rbind,lapply(names(X), function(k) {
    K <- X[[k]]
    do.call(rbind,lapply(names(K), function(w) {
      W <- K[[w]]
      W %>% 
        select(-sum_fitted_p,-MinAge,-min_sum_refitted_p,-diffToMin) %>% 
        rename(sum_fitted_p="sum_refitted_p") %>%
        mutate(window_size=w,
               cluster=k,
               sex=sx) %>%
        select(sex,cluster,window_size,everything())
    }))
  }))
}))
savesource("Fig4e_S6c")
##############################################################################################################################
# @fig3
FigS4a <- da.glm_agebysex.rna %>% 
  mutate(Contrast_label=factor(sub("_",", ",sub("x"," vs ",sub("Age1","HY",sub("Age2","HM",sub("Age3","HO",Contrast_label))))))) %>%
  arrange(-PValue) %>%
  select(-Contrast_label,-Contrast_color,-nhits)
savesource("FigS4a")
##############################################################################################################################
FigS4b <- da.glm_venn_bygene
# R list, saved as compressed R file
save(FigS4b,file = "source_data_final/Fig_S4a.RData")
##############################################################################################################################
# @fig3
FigS4c <- da.glm_agebysex.rna %>% 
  mutate(SexContrast=sub("\\_.*","",Contrast),
         direction=ifelse(logFC<0,"Closing peaks","Opening peaks")) %>%
  inner_join(anno.rna.pbmc %>% select(EnsemblID,GeneName), by="EnsemblID") %>%
  inner_join(with(gsenv,merge(geneset.names.scrnaseq_pbmc_simple_specific,geneset.genes.scrnaseq_pbmc_simple_specific,by="Module.ID")) %>% 
               select(GeneName,Specificity="Module.Name") %>%
               filter(!is.na(Specificity)),
             by="GeneName") %>% 
  mutate(Specificity=ifelse(grepl("^B",Specificity),"CD19",
                            ifelse(grepl("^NK",Specificity),"CD56",
                                   ifelse(grepl("^Mono",Specificity),"CD14",Specificity))),
         Specificity=factor(ifelse(Specificity %in% c(focus.celltypes,"Naive_Tcells","acCD8_Tcells","Tcells"),Specificity,"Other")),
         direction=sub(" in$","",direction)) %>%
  group_by(SexContrast,Specificity) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  mutate(Specificity=factor(Specificity,levels = c(focus.celltypes,"Naive_Tcells","acCD8_Tcells","Tcells","Other")),
         label=factor(paste0(Specificity,"\nn = ",n))) %>%
  mutate(label=factor(label,levels = levels(label)[c(2,4,3,5,7,1,6)])) %>%
  select(sex="SexContrast",direction,Specificity,GeneName,EnsemblID,logFC,PValue,FDR) %>%
  mutate(direction=ifelse(grepl("Closing",direction),"Down","Up"))
savesource("FigS4c")
##############################################################################################################################
# @fig3
FigS4d <- combo.tested.modules.results[[2]] %>%
  filter(Module.Name %in% c(focus.celltypes,"Monocytes","NK_cells","Bcells","acCD8_Tcells","Tcells","Naive_Tcells","Cytotoxic cells","B cells","Inflammation I","Inflammation II","T cells","Myeloid lineage 1","Myeloid lineage 2","B_naive")) %>%
  filter(grepl("Age3xAge1",Contrast_label)) %>%
  mutate(Contrast_color=ifelse(signed.hypergeom.p<0,ifelse(grepl("Age1",.$Contrast_label),colors_age["HY"],colors_age["HM"]),ifelse(grepl("Age3",.$Contrast_label),colors_age["HO"],colors_age["HM"])),
         Contrast_alpha=ifelse(hypergeom.p>combo.threshP,1,0.5),
         Contrast_label=factor(sub("\nAge3xAge1","",Contrast_label)),
         Module.Label=sub(" genes","",as.character(Module.Label))) %>%
  group_by(Module.Label) %>%
  mutate(m=mean(signed.hypergeom.p)) %>%
  ungroup() %>%
  arrange(m) %>%
  mutate(Module.Label=factor(Module.Label,levels = unique(as.character(.$Module.Label)))) %>% 
  left_join(data.frame(.) %>%
              group_by(Module.Label) %>%
              summarize(m=mean(m)) %>%
              mutate(indx=row_number()) %>%
              select(Module.Label,indx),
            by="Module.Label") %>% 
  ungroup() %>% 
  filter(indx<=ntop.gs | indx>=max(indx)-ntop.gs) %>% 
  mutate(signed.hypergeom.p = ifelse(hypergeom.p>maxlogP,sign(signed.hypergeom.p)*maxlogP,signed.hypergeom.p)) %>%
  filter(data.source=="RNA") %>%
  select(celltype="Module.Name",direction="gene.set",gene.count,Total.Freq,hypergeom.p,hypergeom.fdr,signed.hypergeom.p)
# Note: small discrepancies with the data in the figure (number of hits on bars), Figure in ms. must be based on an older version. Pvalues don't
# change (at least perceptibly), so it should be enough to edit the numbers on top of the bars.
savesource("FigS4d")
##############################################################################################################################
# @fig3
gs="Bcells"
transcription <- da.glm_agebysex.rna %>% filter(FDR<=fdr.thresh.rna) %>% inner_join(anno.rna.pbmc,by="EnsemblID")
commonGenes = intersect(gsets_gs[[gs]]$GeneName,transcription$GeneName)
topGSdata <- transcription %>% 
  filter(GeneName %in% commonGenes) %>%
  arrange(GeneName,PValue) %>%
  filter(!duplicated(GeneName)) %>%
  merge(adj.data.rna %>% select(EnsemblID,rnaseq$samps),by="EnsemblID")
FigS4e <- do.call(rbind,lapply(list(Females="F",Males="M"), function(sx) {
  dt <- topGSdata %>% 
    select(GeneName,rnaseq$samps[rnaseq$sex==sx]) %>% 
    column_to_rownames("GeneName") %>% 
    t() %>%
    data.frame() %>%
    rownames_to_column("sampid") %>%
    inner_join(data.frame(age.group=rnaseq$age.group,sex=rnaseq$sex) %>% rownames_to_column("sampid"),by="sampid") %>%
    select(sampid,sex,age.group,everything())
})) %>%
  data.frame(.,row.names = NULL)
# Note: replace sample names with anonymized names
savesource("FigS4e")
##############################################################################################################################
# @fig3
gs="Cytotoxic cells"
transcription <- da.glm_agebysex.rna %>% filter(FDR<=fdr.thresh.rna) %>% inner_join(anno.rna.pbmc,by="EnsemblID")
commonGenes = intersect(gsets_gs[[gs]]$GeneName,transcription$GeneName)
topGSdata <- transcription %>% 
  filter(GeneName %in% commonGenes) %>%
  arrange(GeneName,PValue) %>%
  filter(!duplicated(GeneName)) %>%
  merge(adj.data.rna %>% select(EnsemblID,rnaseq$samps),by="EnsemblID")
FigS4f <- do.call(rbind,lapply(list(Females="F",Males="M"), function(sx) {
  dt <- topGSdata %>% 
    select(GeneName,rnaseq$samps[rnaseq$sex==sx]) %>% 
    column_to_rownames("GeneName") %>% 
    t() %>%
    data.frame() %>%
    rownames_to_column("sampid") %>%
    inner_join(data.frame(age.group=rnaseq$age.group,sex=rnaseq$sex) %>% rownames_to_column("sampid"),by="sampid") %>%
    select(sampid,sex,age.group,everything())
})) %>%
  data.frame(.,row.names = NULL)
# Note: replace sample names with anonymized names
savesource("FigS4f")
##############################################################################################################################
# @fig3
FigS4g <- arc.da.glm_agebysex.nodups %>% 
  filter(abs(DistancetoTSS)<=10e3,ageContrast=="Age3xAge1") %>%
  select(grep("^sex$|^logFC|^PValue|FDR|chr$|start|end|GeneName|Ensembl|Annotation|Distance",colnames(.))) %>%
  select(sex,everything())
savesource("FigS4g")
##############################################################################################################################


# Fig 5 (load Figure4.Rmd Figures5&6.Rmd files)
##############################################################################################################################
# @ fig5&6
Fig5a <- predicted_bysex$common
savesource("Fig5a")
##############################################################################################################################
# @ fig5&6
focus_cs <-c("Monocytes","CD14","NK_cells","CD56","CD4_naive","CD4_memory","CD8_naive","CD8_memory","CD3","B_naive","B_memory","Plasma_cells","CD19")
focus_cs_names <-c("Monocytes","Monocytes","NK cells","NK cells","CD4 naive","CD4 memory","CD8 naive","CD8 memory","T cells","B naive","B memory","Plasma cells","B cells")
Fig5b <- chrspecific.enrichment_bysex.atac$common$fine$all %>%
  select(sex,everything(),-scope) %>%
  filter(celltype %in% focus_cs) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_cs,to = focus_cs_names))
savesource("Fig5b")
##############################################################################################################################
# @ fig5&6
focus_ct <- c("Monocytes","NK_cells","Naive_Tcells","CD4_naive","CD8_naive","Tcells","acCD8_Tcells","Bcells","B_naive","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","CD4 naive","CD8 naive","T cells","Act CD8 T cells","B cells","B cells","Plasma cells")
FigS5c <- geneset.enrichment_bysex.rna$common$dice_major$all %>%
  select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set") %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))
# Note: in ms, an empty row for monocytes is included. This is not in the source data (it is empty)
savesource("FigS5c")
##############################################################################################################################
# @ fig5&6
Fig5d_S7d <- clusters_pooled %>%
  filter(MotifName %in% names(expressed.tf_families))
# Note: this table contains enrichment test results for all tested and expressed TF families. Tables in the ms. have been shortened
# to include only TF families that are enriched in at least one of the contexts in the table (sex/cluster) at FDR<=0.00001
savesource("Fig5d_S7d")
##############################################################################################################################
# @ fig5&6
Fig5e_5f <- discriminantTfs %>% 
  arrange(Cluster,GeneName) %>% 
  filter((xorder<=10 & Cluster==1) | (xorder<=7 & Cluster==3)) %>%
  mutate(Cluster=paste("Cluster ",Cluster)) %>%
  arrange(Cluster,xorder) %>%
  mutate(GeneName=factor(GeneName,levels=unique(.$GeneName))) %>%
  select(sampid,sex,age.group,GeneName,EnsemblID,Cluster,normalized_expression="norm.rna",MotifFamily)
# Replace sample names for anonymized versions
savesource("Fig5d_S7d")
##############################################################################################################################
# @fig4
ctset = "fine"
assign(paste0("logfc_bysex_",ctset),logfc_bysex_byCS[[ctset]]$logfc_bysex %>%
         filter(Specificity.Match %in% cs.included[[ctset]]$Specificity.Alias) %>%
         mutate(Specificity=factor(Specificity,levels = c(cs.included[[ctset]]$Specificity,"Other")),
                Specificity.Match=factor(Specificity.Match,levels = cs.included[[ctset]]$Specificity.Alias)))
FigS5a <- logfc_bysex_fine %>%
  filter(Specificity %in% c("CD4_naive","CD4_memory","CD8_naive","CD8_memory")) %>%
  droplevels() %>%
  select(chr,start,end,GeneName,Annotation,DistancetoTSS,Specificity,Females.atac,Males.atac)
savesource("FigS5a")
##############################################################################################################################
# @fig4
gs = "vp2008"
assign(paste0("logfc_bysex_",gs),logfc_bysex_byGS[[gs]]$logfc_bysex %>%
         left_join(gs.modules.included[[gs]],by="Module.Name") %>%
         filter(Module.Match %in% gs.modules.included[[gs]]$Module.Name) %>%
         mutate(Module.Name=factor(Module.Name,levels = c(gs.modules.included[[gs]]$Module.Name,"Other")),
                Module.Match=factor(Module.Match, levels = gs.modules.included[[gs]]$Module.Name)))
FigS5b <- logfc_bysex_vp2008 %>%
  filter(!Module.Name %in% c("Plasma cells","Myeloid lineage 1","Other")) %>%
  droplevels() %>%
  select(GeneName,Module.Name,Females.rna,Males.rna)
savesource("FigS5b")
##############################################################################################################################
# @fig4
gs = "vp2008"
assign(paste0("logfc_bysex_",gs),logfc_bysex_byGS[[gs]]$logfc_bysex %>%
         left_join(gs.modules.included[[gs]],by="Module.Name") %>%
         filter(Module.Match %in% gs.modules.included[[gs]]$Module.Name) %>%
         mutate(Module.Name=factor(Module.Name,levels = c(gs.modules.included[[gs]]$Module.Name,"Other")),
                Module.Match=factor(Module.Match, levels = gs.modules.included[[gs]]$Module.Name)))
FigS5c <- logfc_bysex_vp2008 %>%
  filter(!Module.Name %in% c("Plasma cells","Myeloid lineage 1","Other")) %>%
  droplevels() %>%
  select(chr,start,end,GeneName,Annotation,DistancetoTSS,Module.Name,Females.atac,Males.atac)
savesource("FigS5c")
##############################################################################################################################


# Fig 6 (load Figures5&6.Rmd Figure7.Rmd files)
##############################################################################################################################
#@ fig7
Fig6a <- da.glm_sexbyage.atac %>% 
  mutate(Contrast_label=plyr::mapvalues(Contrast_label,from=levels(.$Contrast_label), to = sub("_",", ",sub("x"," vs ",sub("Age1","HY",sub("Age2","HM",sub("Age3","HO",levels(.$Contrast_label)))))))) %>%
  select(-Contrast_color,-nhits,-Contrast_label) %>%
  mutate(Age=sub("\\_MxF","",Contrast)) %>%
  select(Age,everything(),-Contrast)
savesource("Fig6a")
##############################################################################################################################
#@ fig7
Fig6b <- da.glm_sexbyage.rna %>% 
  mutate(Contrast_label=plyr::mapvalues(Contrast_label,from=levels(.$Contrast_label), to = sub("_",", ",sub("x"," vs ",sub("Age1","HY",sub("Age2","HM",sub("Age3","HO",levels(.$Contrast_label)))))))) %>%
  select(-Contrast_color,-nhits,-Contrast_label) %>%
  mutate(Age=sub("\\_MxF","",Contrast)) %>%
  select(Age,GeneName,everything(),-Contrast)
savesource("Fig6b")
##############################################################################################################################
#@ fig7
Fig6c <- rbind(read.table("../Tables/unused/fig7_atac_enrichment_scrnaseq_pbmc_simple_specific.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F) %>%
                 select(celltype="Module.Name",sex_bias="gene.set",everything(),-dummypos,-Module.Label,-Contrast_label) %>%
                 mutate(source="ATAC-seq"),
               read.table("../Tables/unused/fig7_rna_enrichment_scrnaseq_pbmc_simple_specific.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F) %>%
                 select(celltype="Module.Name",sex_bias="gene.set",everything(),-dummypos,-Module.Label,-Contrast_label) %>%
                 mutate(source="RNA-seq"))
savesource("Fig6c")
##############################################################################################################################
# @ fig5&6
focus_ct <-c("Monocytes","NK_cells","Naive_Tcells","Tcells","acCD8_Tcells","Bcells","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","T cells","Act CD8 T cells","B cells","Plasma cells")
FigS6a1 <- rbind(geneset.enrichment_bysex.atac$females$scrnaseq_pbmc_simple_specific$test$all %>%
                   select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set"),
                 geneset.enrichment_bysex.atac$males$scrnaseq_pbmc_simple_specific$test$all %>%
                   select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set")) %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))

focus_cs <-c("Monocytes","CD14","NK_cells","CD56","CD4_naive","CD4_memory","CD8_naive","CD8_memory","CD3","B_naive","B_memory","Plasma_cells","CD19")
focus_cs_names <-c("Monocytes","Monocytes","NK cells","NK cells","CD4 naive","CD4 memory","CD8 naive","CD8 memory","T cells","B naive","B memory","Plasma cells","B cells")
savesource("FigS6a1")
FigS6a2 <- rbind(chrspecific.enrichment_bysex.atac$females$fine$all %>%
                   select(sex,everything(),-scope),
                 chrspecific.enrichment_bysex.atac$males$fine$all %>%
                   select(sex,everything(),-scope)) %>%
  filter(celltype %in% focus_cs) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_cs,to = focus_cs_names))
savesource("FigS6a2")
##############################################################################################################################
# @ fig5&6
focus_ct <- c("Monocytes","NK_cells","Naive_Tcells","CD4_naive","CD8_naive","Tcells","acCD8_Tcells","Bcells","B_naive","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","CD4 naive","CD8 naive","T cells","Act CD8 T cells","B cells","B cells","Plasma cells")
FigS6b <- rbind(geneset.enrichment_bysex.rna$females$dice_major$all %>%
                 select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set"),
               geneset.enrichment_bysex.rna$males$dice_major$all %>%
                 select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set")) %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))
savesource("FigS6b")
##############################################################################################################################


# Fig 7 (load Figures5&6.Rmd)
##############################################################################################################################
# @ fig5&6
FigS7a <- predicted_bysex.rna$common
savesource("FigS7a")
##############################################################################################################################
focus_ct <-c("Monocytes","NK_cells","Naive_Tcells","Tcells","acCD8_Tcells","Bcells","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","T cells","Act CD8 T cells","B cells","Plasma cells")
FigS7b <- geneset.enrichment_bysex.atac$common$scrnaseq_pbmc_simple_specific$test$all %>%
  select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set") %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))
savesource("FigS7b")
##############################################################################################################################
focus_ct <-c("Monocytes","NK_cells","Naive_Tcells","CD4_naive","CD8_naive","Tcells","acCD8_Tcells","Bcells","B_naive","Plasma_cells")
focus_ct_names <-c("Monocytes","NK cells","Naive T cells","CD4 naive","CD8 naive","T cells","Act CD8 T cells","B cells","B cells","Plasma cells")
FigS7c <- geneset.enrichment_bysex.rna$common$scrnaseq_pbmc_simple_specific$all %>%
  select(sex,everything(),-Module.Label,-schema,-dummypos,cluster="gene.set") %>%
  rename(celltype="Module.Name") %>%
  filter(celltype %in% focus_ct) %>%
  mutate(celltype=plyr::mapvalues(celltype,from = focus_ct,to = focus_ct_names))
# Note: in ms, 4 empty rows are included. These are not in the source data (they're empty)
savesource("FigS7c")
##############################################################################################################################
# @ fig5&6
FigS7dr <- adjusted.top_bysex$common %>% # rna chunk, ntop = 8
  select(GeneName,sampid,sex,Cluster,age.group,age,Specificity,adjusted_expression="adj.score")
savesource("FigS7dr")
FigS7da <- adjusted.top_bysex$common %>% # atac chunk, use fine-grained states, ntop = 12
  select(peakco,sampid,sex,Cluster,age.group,age,Specificity,adjusted_accessibility="adj.score")
savesource("FigS7da")
##############################################################################################################################

# Fig 8 (load Figure7.Rmd)
##############################################################################################################################
FigS8a <- rbind(read.table("../Tables/unused/fig7_atac_enrichment_vp2008.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F) %>%
                  select(celltype="Module.Name",sex_bias="gene.set",everything(),-dummypos,-Module.Label,-Contrast_label) %>%
                  mutate(source="ATAC-seq"),
                read.table("../Tables/unused/fig7_rna_enrichment_vp2008.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F) %>%
                  select(celltype="Module.Name",sex_bias="gene.set",everything(),-dummypos,-Module.Label,-Contrast_label) %>%
                  mutate(source="RNA-seq"))
# Consider renaming modules to match names in other figures
savesource("FigS8a")
##############################################################################################################################
FigS8b1 <- da.chromhmm.spec_sexbyage %>%
  filter(FDR<=fdr.thresh) %>%
  filter(AgeContrast=="HO") %>%
  select(peakco,AgeContrast,direction,Specificity)
savesource("FigS8b1")
FigS8b2 <- da.scrnaseq.spec_sexbyage %>%
  filter(FDR<=fdr.thresh.rna) %>%
  filter(AgeContrast=="HO") %>%
  select(GeneName,EnsemblID,AgeContrast,direction,Specificity)
savesource("FigS8b2")
##############################################################################################################################
FigS8c1 <- ca.chromhmm.spec_sexbyage %>%
  group_by(AgeContrast,Specificity) %>%
  mutate(n=n()) %>%
  group_by(AgeContrast,Specificity,direction) %>%
  mutate(ndir=n(),
         ndir=ifelse(!duplicated(ndir),ndir,NA)) %>%
  ungroup() %>%
  filter(AgeContrast=="HO") %>%
  droplevels() %>%
  mutate(chromHMMstate="Cell-specific\nregulatory elements",
         Specificity=factor(sub("Plasma","B_Plasma",Specificity),levels = rev(sub("Plasma","B_Plasma",focus_cs)))) %>%
  left_join(data.frame(.) %>% 
              dcast(.,peakco+AgeContrast+Specificity~sex,value.var="mean.score") %>%
              group_by(AgeContrast,Specificity) %>%
              summarize(WilcoxP=wilcox.test(Females,Males)$p.value) %>%
              ungroup(),
            by=c("AgeContrast","Specificity")) %>%
  mutate(WilcoxP_label=ifelse(WilcoxP>0.01,"n.s.",ifelse(WilcoxP<0.0001,"***",ifelse(WilcoxP<0.001,"**","*"))),
         label=paste0(Specificity,WilcoxP_label,", n = ",n)) %>%
  select(peakco,AgeContrast,sex,mean_normalized_accessibility="mean.score",logFC,logCPM,LR,PValue,FDR,Specificity)
savesource("FigS8c1")
FigS8c2 <- ex.scrnaseq.spec_sexbyage <- data.frame(rnaseq$da_byage$Age3_MxF$da.list$adj) %>% 
  rownames_to_column("EnsemblID") %>%
  melt(.,id.vars="EnsemblID",variable.name="sampid",value.name="adj.score") %>%
  mutate(sampid=as.character(sampid)) %>%
  inner_join(with(rnaseq,data.frame(sampid,sex,age.group,age,stringsAsFactors=F)),by="sampid") %>%
  group_by(EnsemblID) %>%
  mutate(zscore=znorm(adj.score)) %>%
  group_by(EnsemblID,age.group,sex) %>%
  summarize(mean.score=mean(zscore)) %>%
  ungroup() %>%
  mutate(sex=plyr::mapvalues(sex,from = c("F","M"),to = c("Females","Males"))) %>%
  dplyr::rename(AgeContrast="age.group") %>%
  inner_join(da.scrnaseq.spec_sexbyage,by=c("EnsemblID","AgeContrast")) %>%
  filter(!is.na(Specificity)) %>%
  mutate(Specificity=ifelse(Specificity %in% c("DCs","pDCs"),"Dendritic_cells",Specificity),
         Specificity=ifelse(Specificity %in% focus_sc,Specificity,"Other")) %>%
  select(GeneName,EnsemblID,AgeContrast,sex,mean_normalized_expression="mean.score",logFC,logCPM,LR,PValue,FDR,Specificity)
savesource("FigS8c2")
##############################################################################################################################


savesource <- function(x) write.table(get(x),file = paste0("source_data_final/",sub("Fig","Fig_",x),".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
