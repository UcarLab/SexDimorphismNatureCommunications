# Preps (loads, filters, annotates) ATACseq base data for downstream analyses

library(edgeR)
library(reshape2)
library(biomaRt)
library(GenomicRanges)
library(dplyr)

rm(list = ls())

andred <- function(x) {Reduce("&",x)}
# Peak blacklist for hg19, obtained by merging DNAse-seq blacklisted regions downloaded from UCSC genome browser
blacklist <- GRanges(read.table("../../Datasets/genome/hg19.blacklists_merged.bed",header = F,sep = "\t",quote = "",col.names = c("chr","start","end","annotation")))
  
statsfile = "filtered_global.stats.RData"
load(statsfile) # metadata variables for samples used in the study, as separate vectors: samps,sampid,  age.group, sex, age, age.group, season, batch, batchdate, nreps, frailty, etc

# narrowPeaks data
# Here, pooled_rawcounts.txt contains a matrix formed by pasting raw read counts of all samples at the same consensus peaks, with samples as columns and peaks as rows
atac.pbmc.narrow <- read.table("data/narrowPeaks/pooled_rawcounts.txt",header = T,sep = "\t",quote = "") %>%
    arrange(chr,start,end) %>%
    setNames(.,sub("_PBMC","",colnames(.)))
atac.narrow.gr <- GRanges(atac.pbmc.narrow)
atac.narrow.gr.whitelisted.idx <- setdiff(seq_len(length(atac.narrow.gr)),findOverlaps(atac.narrow.gr,blacklist,ignore.strand=T)@from)
atac.narrow.gr.whitelisted <- atac.narrow.gr[atac.narrow.gr.whitelisted.idx]

atac.basepeaks.pbmc <- data.frame(atac.narrow.gr.whitelisted,stringsAsFactors = F) %>% 
  rename(chr="seqnames") %>%
  mutate(peakco=paste(chr,start,end,sep = "_"),chr=as.character(chr)) %>%
  data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F) %>%
  select(chr,start,end,peakco,nhits,sampid)
write.table(atac.basepeaks.pbmc,"./data/narrowPeaks/narrow_raw_whitelisted.txt",quote = F,sep = "\t",row.names = F,col.names = T)

libsize.atac.pbmc <- colSums(atac.basepeaks.pbmc[,sampid])

# Reads and processes peakStats files
peaksignals.pbmc <- (read.table("./data/narrowPeaks/peakStats/narrow_peakstats_signals.txt",header = T,sep = "\t",stringsAsFactors = F,row.names = "peakco") %>%
  setNames(.,sub("signal_","",colnames(.))) %>%
  select(chr,start,end,sampid))[atac.basepeaks.pbmc$peakco,]
peakcalls.pbmc <- cbind(peaksignals.pbmc %>% select(chr,start,end),1*(peaksignals.pbmc[,sampid]>0)) %>%
  mutate(HY.counts=rowSums(.[sampid[age.group=="HY"]]),
         HM.counts=rowSums(.[sampid[age.group=="HM"]]),
         HO.counts=rowSums(.[sampid[age.group=="HO"]]),
         Total.counts=rowSums(.[sampid]))
peakqvalues.pbmc <- (read.table("./data/narrowPeaks/peakStats/narrow_peakstats_qvalues.txt",header = T,sep = "\t",stringsAsFactors = F,row.names = "peakco") %>%
  setNames(.,sub("qvalue_","",colnames(.))) %>%
  select(chr,start,end,sampid))[atac.basepeaks.pbmc$peakco,] %>%
  mutate(maxq=apply(.[,sampid],1,max)) %>%
  data.frame(.,row.names = atac.basepeaks.pbmc$peakco,stringsAsFactors = F)
atac.allpeaks.pbmc <- atac.basepeaks.pbmc %>%
  filter(peakcalls.pbmc$Total.counts>=2,
         peakqvalues.pbmc$maxq>-log10(0.01),
         !chr %in% c("chrX","chrY","chrM"),
         apply(.[,sampid],1,max)>20, 
         apply(cpm(.[,sampid],lib.size = libsize.atac.pbmc),1,max)<500)
write.table(atac.allpeaks.pbmc,"./data/narrowPeaks/narrow_raw_whitelisted_filtered.txt",row.names = F, col.names = T, sep = "\t", quote = F)

### Homer and chromHMM annotations
# To get this file run narrow_raw_whitelisted_filtered.txt, saved above, through add_homer_annotations.sh
atac.allpeaks.anno <- read.table("./data/narrowPeaks/narrow_raw_whitelisted_filtered_annotated.txt", header = T, sep = "\t", quote = "",stringsAsFactors = F) %>%
  mutate(peakco=paste(chr,start,end,sep = "_")) %>%
  data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F)

# To get these files run narrow_raw_whitelisted_filtered.txt, saved above, through add_chromhmm_annotations.sh
chrmhmm.sets <- setNames(as.list(unique(sort(sub(".*filtered\\_","",sub("\\_chromHMM.txt","",dir(path = "./data/narrowPeaks/chromHMM/")[grepl("\\_chromHMM\\.txt",dir(path = "./data/narrowPeaks/chromHMM/"))]))))),
                         c("B_memory",
                           "CD14",
                           "CD15",
                           "CD19",
                           "CD34",
                           "CD4_memory",
                           "CD4_naive",
                           "Th",
                           "Th17",
                           "Treg",
                           "CD56",
                           "CD8_naive",
                           "CD8_memory",
                           "CD34_mobilized",
                           "PBMC",
                           "B_germinal_center",
                           "B_naive",
                           "B_memory_nonclass_switched",
                           "Plasma_cells",
                           "Macrophages_inflammatory",
                           "Erythroblasts",
                           "Megakaryocytes",
                           "Macrophages_alt_activated"))[c("CD14",
                                                           "CD19",
                                                           "B_naive",
                                                           "B_memory",
                                                           "Plasma_cells",
                                                           "CD4_naive",
                                                           "CD4_memory",
                                                           "Th",
                                                           "Th17",
                                                           "Treg",
                                                           "CD56",
                                                           "CD8_naive",
                                                           "CD8_memory",
                                                           "Erythroblasts",
                                                           "Megakaryocytes",
                                                           "PBMC")]

## ChromHMM + priority selection to retain a unique chromHMM state per peak
proximalThresh = 1000 # peaks within this distance (in bp) from TSS are considered "proximal" for resolving Tss vs Enh conflicting anntoations: 
                      # a peak annotated as both is considered a Tss if proximal and Enh if distal (i.e. not proximal) by comparison to this threshold
chmm.prior <- data.frame(
  Priority=1:30,
  chromHMMstate.proximal=c("TssA","TssA1","TssA2","TssWk","EnhA1","EnhA2","EnhG1","EnhG2","Enh","EnhG","TssBiv","EnhWk","EnhBiv","ReprPC","ReprPCWk","Repr1","Repr2",
                           "TssFlnk","TssFlnkU","TssFlnkD","TssAFlnk","TssADistal","BivTssFlnk","Tx","TxElong","TxWk","TxFlnk","ZNF_Rpts","Het","Quies"),
  chromHMMstate.distal=c("EnhA1","EnhA2","EnhG1","EnhG2","Enh","EnhG","EnhWk","EnhBiv","TssA","TssA1","TssA2","TssWk","TssBiv","ReprPC","ReprPCWk","Repr1","Repr2",
                         "TssFlnk","TssFlnkU","TssFlnkD","TssAFlnk","TssADistal","BivTssFlnk","Tx","TxElong","TxWk","TxFlnk","ZNF_Rpts","Het","Quies"),
  stringsAsFactors = F
)

peaks.hmm <- lapply(chrmhmm.sets, function(cs) {
  hmm <- read.table(paste("./data/narrowPeaks/chromHMM/narrow_raw_whitelisted_filtered",cs,"chromHMM.txt",sep = "_"), header = T, sep = "\t", quote = "",stringsAsFactors = F) %>%
    mutate(peakco=paste(chr,start,end,sep = "_")) %>%
    left_join(atac.allpeaks.anno %>% select(peakco,DistancetoTSS),by="peakco") %>%
    mutate(isProximal=ifelse(abs(DistancetoTSS<=proximalThresh),TRUE,FALSE))
  hmmPeaks <- (rbind(hmm %>% 
                       filter(isProximal) %>% 
                       left_join(chmm.prior %>% select(chromHMMstate.proximal,Priority) %>% rename(chromHMMstate="chromHMMstate.proximal"),by = "chromHMMstate"),
                     hmm %>% 
                       filter(!isProximal) %>% 
                       left_join(chmm.prior %>% select(chromHMMstate.distal,Priority) %>% rename(chromHMMstate="chromHMMstate.distal"),by = "chromHMMstate")) %>%
                 arrange(peakco,Priority) %>%
                 filter(!duplicated(peakco)) %>%
                 select(chr,start,end,peakco,chromHMMstate) %>%
                 data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F))[atac.allpeaks.pbmc$peakco,] %>%
    mutate(chromHMMsimple=sub("^Het.*","Others",
                              sub("^Repr[1,2]","Others",
                                  sub("^ZNF.*","Others",
                                      sub("^Tx.*","Tx",
                                          sub(".*Tss.*","Tss",
                                              sub("ReprPC.*","ReprPC",
                                                  sub(".*Enh.*","Enh",
                                                      chromHMMstate)))))))) %>%
    rename_at(vars(grep("chromHMM",colnames(.))),funs(sub("$",paste0("_",names(chrmhmm.sets[chrmhmm.sets==cs])),.)))
  write.table(hmmPeaks,paste("./data/narrowPeaks/chromHMM/narrow_whitelisted_filtered",cs,"chromHMM_single.txt",sep = "_"),sep = "\t",col.names = T, row.names = F, quote = F)
  return(hmmPeaks)
})

# Since B cell states were derived from hg38 data and had to be liftover to hg19, not all segments carried forward. This step replaces these NA peaks with CD19 states
for (ct in names(peaks.hmm)[grep("^B_|^Plasma",names(peaks.hmm))]) {
  nna = which(is.na(peaks.hmm[[ct]]$start))
  peaks.hmm[[ct]][nna,] <- peaks.hmm$CD19[nna,]
}

atac.allpeaks.hmm.pbmc <- Reduce(function(x,y) merge(x,y,by=c("chr","start","end","peakco")), peaks.hmm)
write.table(atac.allpeaks.hmm.pbmc,"./data/narrowPeaks/narrow_whitelisted_filtered_chromHMM_single.txt",sep = "\t",col.names = T, row.names = F, quote = F)

# Finding cell type-specific _active_ enhancers and promoters
act_Enh <- c("EnhA1","EnhA2","EnhA","EnhG1","EnhG2","Enh")
act_Tss <- c("TssA","TssA1","TssA2")
act_RE <- c(act_Enh,act_Tss)

# The following pipe establishes cell type specificity of each consensus peak at various scales, described in Supplementary Materials
atac.allpeaks.cspec_hmm.pbmc <- atac.allpeaks.hmm.pbmc %>% 
  select(grep("peakco|chromHMMstate",colnames(.))) %>%
  select(-chromHMMstate_PBMC,-chromHMMstate_CD19) %>%
  rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
  mutate(EnhA_specificity_granular_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
         TssA_specificity_granular_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
         RE_specificity_granular_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
  left_join(data.frame(.) %>% 
              filter(EnhA_specificity_granular_n==1) %>%
              mutate(EnhA_specificity_granular=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
              select(peakco,EnhA_specificity_granular),
            by="peakco") %>%
  left_join(data.frame(.) %>% 
              filter(TssA_specificity_granular_n==1) %>%
              mutate(TssA_specificity_granular=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
              select(peakco,TssA_specificity_granular),
            by="peakco") %>%
  left_join(data.frame(.) %>% 
              filter(RE_specificity_granular_n==1) %>%
              mutate(RE_specificity_granular=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
              select(peakco,RE_specificity_granular),
            by="peakco") %>%
  select(peakco,EnhA_specificity_granular,TssA_specificity_granular,RE_specificity_granular) %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC,-chromHMMstate_Th,-chromHMMstate_Th17,-chromHMMstate_Treg) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              select(grep("CD19",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_fine_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_fine_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_fine_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_fine_n==1) %>%
                          mutate(EnhA_specificity_fine=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_fine),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_fine_n==1) %>%
                          mutate(TssA_specificity_fine=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_fine),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_fine_n==1) %>%
                          mutate(RE_specificity_fine=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_fine),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_fine,TssA_specificity_fine,RE_specificity_fine),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(CD19=data.frame(.) %>% select(grep("^B|CD19|Plasma",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     CD4T=data.frame(.) %>% select(grep("CD4|^Th|^Treg",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     CD8T=data.frame(.) %>% select(grep("CD8",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     CD14=data.frame(.) %>% select(grep("CD14",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other")))) %>%
              select(grep("CD4\\_|CD8\\_|B\\_|Plasma|^Th|^Treg",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_gross_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_gross_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_gross_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_gross_n==1) %>%
                          mutate(EnhA_specificity_gross=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_gross),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_gross_n==1) %>%
                          mutate(TssA_specificity_gross=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_gross),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_gross_n==1) %>%
                          mutate(RE_specificity_gross=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_gross),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_gross,TssA_specificity_gross,RE_specificity_gross),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC,-chromHMMstate_Th,-chromHMMstate_Th17,-chromHMMstate_Treg) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(naive_B=data.frame(.) %>% select(grep("B_naive",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     memory_B=data.frame(.) %>% select(grep("B_memory",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     naive_T=data.frame(.) %>% select(grep("CD.\\_naive",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     memory_T=data.frame(.) %>% select(grep("CD.\\_memory",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     Monocytes=data.frame(.) %>% select(grep("CD14",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other")))) %>%
              select(grep("\\_naive|\\_memory|^CD1",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_state_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_state_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_state_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_state_n==1) %>%
                          mutate(EnhA_specificity_state=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_state),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_state_n==1) %>%
                          mutate(TssA_specificity_state=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_state),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_state_n==1) %>%
                          mutate(RE_specificity_state=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_state),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_state,TssA_specificity_state,RE_specificity_state),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(CD3=data.frame(.) %>% select(grep("CD4|CD8|^Th|^Treg",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     CD19=data.frame(.) %>% select(grep("CD19|^B|Plasma",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     CD14=data.frame(.) %>% select(grep("CD14",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other")))) %>%
              select(grep("CD4|CD8|^B|Plasma|^Th|^Treg",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_sublineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_sublineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_sublineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_sublineage_n==1) %>%
                          mutate(EnhA_specificity_sublineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_sublineage),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_sublineage_n==1) %>%
                          mutate(TssA_specificity_sublineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_sublineage),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_sublineage_n==1) %>%
                          mutate(RE_specificity_sublineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_sublineage),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_sublineage,TssA_specificity_sublineage,RE_specificity_sublineage),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(Myeloid=data.frame(.) %>% select(grep("CD14|Erythro|Megakaryo",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     Lymphoid=data.frame(.) %>% select(grep("CD4|CD8|^B|Plasma|CD56|^Th|^Treg",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other")))) %>%
              select(grep("CD4|CD8|^B|Plasma|^CD1|CD56|Erythro|Megakaryo|^Th|^Treg",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_lineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_lineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_lineage_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_lineage_n==1) %>%
                          mutate(EnhA_specificity_lineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_lineage),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_lineage_n==1) %>%
                          mutate(TssA_specificity_lineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_lineage),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_lineage_n==1) %>%
                          mutate(RE_specificity_lineage=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_lineage),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_lineage,TssA_specificity_lineage,RE_specificity_lineage),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(Innate=data.frame(.) %>% select(grep("CD14|CD56",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other"))),
                     Adaptive=data.frame(.) %>% select(grep("CD4|CD8|^B|Plasma|^Th|^Treg",colnames(.))) %>% apply(.,1,function(x) ifelse(any(x %in% act_Enh),"EnhA",ifelse(any(x %in% act_Tss),"TssA","Other")))) %>%
              select(grep("CD4|CD8|^B|Plasma|^CD1|CD56|^Th|^Treg",colnames(.),invert=T)) %>%
              mutate(EnhA_specificity_function_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_function_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_function_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_function_n==1) %>%
                          mutate(EnhA_specificity_function=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_function),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_function_n==1) %>%
                          mutate(TssA_specificity_function=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_function),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_function_n==1) %>%
                          mutate(RE_specificity_function=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_function),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_function,TssA_specificity_function,RE_specificity_function),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(peakco,chromHMMstate_CD4_naive,chromHMMstate_Th,chromHMMstate_Th17,chromHMMstate_Treg) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(EnhA_specificity_CD4_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Enh)),
                     TssA_specificity_CD4_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_Tss)),
                     RE_specificity_CD4_n=data.frame(.) %>% apply(.,1,function(x) sum(x %in% act_RE))) %>%
              left_join(data.frame(.) %>% 
                          filter(EnhA_specificity_CD4_n==1) %>%
                          mutate(EnhA_specificity_CD4=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Enh)])) %>%
                          select(peakco,EnhA_specificity_CD4),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(TssA_specificity_CD4_n==1) %>%
                          mutate(TssA_specificity_CD4=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_Tss)])) %>%
                          select(peakco,TssA_specificity_CD4),
                        by="peakco") %>%
              left_join(data.frame(.) %>% 
                          filter(RE_specificity_CD4_n==1) %>%
                          mutate(RE_specificity_CD4=data.frame(.) %>% select(-peakco) %>% apply(.,1,function(x) colnames(.)[which(x %in% act_RE)])) %>%
                          select(peakco,RE_specificity_CD4),
                        by="peakco") %>%
              select(peakco,EnhA_specificity_CD4,TssA_specificity_CD4,RE_specificity_CD4),
            by="peakco") %>%
  left_join(atac.allpeaks.hmm.pbmc %>% 
              select(grep("peakco|chromHMMstate",colnames(.))) %>%
              select(-chromHMMstate_PBMC,-chromHMMstate_Th,-chromHMMstate_Th17,-chromHMMstate_Treg) %>%
              rename_all(funs(sub("chromHMMstate\\_","",.))) %>%
              mutate(
                chromHMM_Diversity=apply(.,1,function(x) {
                  p <- prop.table(table(as.character(x[-1])))
                  H=-sum(p*log(p))
                  return(H)}),
                N=apply(.,1,function(x) length(unique(as.character(x[-1])))),
                chromHMM_Evenness=ifelse(N==1,1,chromHMM_Diversity/log(N)),
                chromHMM_invariant_state=ifelse(chromHMM_Diversity==0,CD14,NA)) %>%
              select(peakco,chromHMM_Diversity,chromHMM_Evenness,chromHMM_invariant_state),
            by="peakco")

# Putting everything together and save annotation frame
atac.allpeaks.anno.pbmc <- atac.allpeaks.anno %>% 
  left_join(atac.allpeaks.hmm.pbmc,by = c("chr","start","end","peakco")) %>%
  left_join(atac.allpeaks.cspec_hmm.pbmc,by = c("peakco")) %>%
  data.frame(.,row.names = paste(.$chr,.$start,.$end,sep = "_"),stringsAsFactors = F)
atac.allpeaks.anno.pbmc <- atac.allpeaks.anno.pbmc[atac.allpeaks.pbmc$peakco,]
write.table(atac.allpeaks.anno.pbmc,"./data/narrowPeaks/narrow_whitelisted_filtered_annotated.txt",sep = "\t",col.names = T, row.names = F, quote = F)

saveobj <- ls()[grep("pbmc",ls())]
save(list = saveobj,file = "./data/narrowPeaks/pbmc_narrowPeaks.RData")

