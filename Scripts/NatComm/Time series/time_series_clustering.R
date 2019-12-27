# Function to perform k-means clustering of time series data estimated using auto-arima, and scripts used
# to apply this function to ATAC-seq data. Same scripts are used for RNAseq data, with minor modifications

Temp=new.env(parent = .GlobalEnv)

library(parallel)
library(reshape2)
library(dplyr)

###################################################################################################################################################################
## Function: KMA analysis of time series trajectories
# Includes a selection of optimal k (n clusters) based on elbow detection of accumulated variance within clusters relative to total variance. 
# To access this, nrand has to be set to > 0, which prompts the calculation of clustering statistics at each k between 1 and max_k.
# If nrand=0, max_k is the only k used.
# recluster_k is used to collapse k-means clusters by joining the most similar clusters. This is done by applying hierarchical clustering to the
# k cluster medians and using cutree() to join nearest clusters.
# To determine "significance" of a k, results based on a sequence of k values are compared to a similar sequence based on randomly selected data
# rows. The highest k for which the variance captured by clusters relative to the variance outside of clusters exceeds the same value computed from
# nrand random permutations of the data is chosen as the optimal k value.
# kmeans_on_zscores determines whether to normalize (z-score) time series data prior to analysis.
###################################################################################################################################################################

kbins_peaks <- function(peaks_ts=NULL,
                        binage=NULL,
                        recluster_k = NULL,
                        max_k=2,
                        nrand=1000,
                        p.rand=0.05,
                        seed=1,
                        kmeans_on_zscores=FALSE,
                        dataType=c("atac","rna"),
                        n.cores=2) {
  
  require(cluster)
  require(MKmisc)
  require(reshape2)
  require(parallel)
  require(dplyr)
  
  if (any(names(binage)!=names(peaks_ts))) {
    cat("Warning: names of BINAGE construct do not match names of PEAKS_TS matrix. Analysis cannot continue.")
    return(NULL)
  }
  if (is.null(peaks_ts)) {
    cat("Warning: input data (peaks_ts) matrix is empty. Analysis cannot continue.")
    return(NULL)
  }
  if (recluster_k>max_k) {
    cat("Warning: recluster k (",recluster_k,") cannot be larger than max k (",max_k,"). It has been reduced to match this value.")
    recluster_k=max_k
  }
  
  random.distpeaks.wss <- function(tsmx,kmx,nrand,n.cores) {
    set.seed = seed
    clust.split <- lapply(kmx,`[[`,"cluster")
    peakarray <- data.frame(peakco=rownames(tsmx))
    tot.wss <- mclapply(seq_len(nrand),function(i) {
      random.clust <- lapply(clust.split,function(x) cbind(peakarray,cluster=sample(x,length(x),replace = F)))
      tot.wss <- lapply(random.clust, function(x) sum(unlist(lapply(seq_len(max(x$cluster)), function(j) kmeans(tsmx[x$cluster==j,],centers = 1)$totss))))
      return(tot.wss)
    },mc.cores = n.cores)
    tot.wss.df <- setNames(do.call(data.frame,lapply(tot.wss,unlist)),seq_len(nrand))
    return(tot.wss.df)
  }
  zscores <- function(x) (x-mean(x))/sd(x)
  
  peaks_ts.nonzero <- data.frame(peaks_ts,check.names = F)
  if (kmeans_on_zscores) {
    # Kmeans computed on zscores
    cat("Computing",ifelse(dataType=="atac","peak","transcript"),"z-scores\n")
    peaks_ts.nonzero <- t(apply(peaks_ts.nonzero,1,function(x) (x-mean(x))/sd(x)))
  }
  # Correlation and distance matrices
  corage <- cor(peaks_ts.nonzero)
  cordistage <- corDist(t(peaks_ts.nonzero))
  distage <- dist(t(peaks_ts.nonzero),method = "euclidean")
  bincuts <- which(!!diff(binage))
  
  if (nrand>0) {
    if (p.rand < 1/nrand) cat("Warning:",paste0("the selected N of random replicates (",nrand,") is too low to see significance at the chosen level (",p.rand,"). A minimum of N=",ceiling(1/p.rand)," is required to see significance at that level.\n"))
    
    # "Observed" clusters
    cat("Computing k-means on peak time series\n")
    distpeaks.kmeans <- mclapply(seq_len(max_k),function(x) kmeans(peaks_ts.nonzero,centers=x,iter.max = 100,nstart = 1000),
                                 mc.cores = n.cores)
    distpeaks.kmeans.wss <- data.frame(tot.wss=unlist(lapply(distpeaks.kmeans,function(x) x$tot.withinss)),row.names = paste("k",1:max_k,sep = ""))
    distpeaks.kmeans.wss.prop <- distpeaks.kmeans.wss/max(distpeaks.kmeans.wss)
    
    # First differentials to estimate elbow location
    delta.distpeaks.kmeans.wss <- data.frame(Obs=abs(diff(distpeaks.kmeans.wss$tot.wss)),row.names = paste("k",seq_len(max_k)[-1],sep = ""))
    delta.distpeaks.kmeans.wss.prop <- data.frame(Obs=abs(diff(distpeaks.kmeans.wss.prop$tot.wss)),row.names = paste("k",seq_len(max_k)[-1],sep = ""))
    rm("distpeaks.kmeans.wss");gc()
    
    # Randomized clusters
    cat("Computing k-means on randomized data (n=",nrand,")\n")
    random.distpeaks.kmeans.wss <- random.distpeaks.wss(peaks_ts.nonzero,distpeaks.kmeans,nrand,n.cores)
    rm("distpeaks.kmeans");gc()
    random.distpeaks.kmeans.wss.prop <- data.frame(random.distpeaks.kmeans.wss %>% mutate_all(funs(./max(.))),row.names = paste0("k",1:max_k),check.names = F)
    
    # First differentials
    delta.rand.distpeaks.kmeans.wss <- data.frame(rbind(vector("numeric",nrand),apply(random.distpeaks.kmeans.wss,2,function(x) abs(diff(as.numeric(x))))),row.names = paste0("k",seq_len(max_k)),check.names = F)
    rm("random.distpeaks.kmeans.wss");gc()
    delta.rand.distpeaks.kmeans.wss.prop <- data.frame(rbind(vector("numeric",nrand),apply(random.distpeaks.kmeans.wss.prop,2,function(x) abs(diff(as.numeric(x))))),row.names = paste0("k",seq_len(max_k)),check.names = F)
    
    # Probability value for each k, returns cluster number estimate
    prob.distpeaks.kmeans.wss <- data.frame(Prob=((1+rowSums(apply(delta.rand.distpeaks.kmeans.wss,2,function(x) x>c(1,delta.distpeaks.kmeans.wss$Obs))))/(nrand+1))) %>%
      slice(-1) %>%
      mutate(Signif=Prob<=p.rand) %>%
      data.frame(row.names = rownames(delta.distpeaks.kmeans.wss))
    K <- sum(prob.distpeaks.kmeans.wss$Signif)+1
    rm(list = c("delta.distpeaks.kmeans.wss","delta.rand.distpeaks.kmeans.wss"));gc()
    
    # Prep of visualization of results of randomization tests
    cat("Preparing output for visualization of randomization test results\n")
    kopt.peaks.wss.prop <- melt(merge(distpeaks.kmeans.wss.prop,random.distpeaks.kmeans.wss.prop,by="row.names"),id.vars = "Row.names",variable.name = "Run",value.name = "Prop") %>%
      mutate(Cluster=factor(Row.names,levels = paste("k",1:max_k,sep = "")),Run=ifelse(Run=="tot.wss","obs","rand")) %>%
      select(-Row.names)
    rm(list = c("random.distpeaks.kmeans.wss.prop","distpeaks.kmeans.wss.prop"));gc()
    delta.kopt.peaks.wss.prop <- merge(melt(merge(delta.distpeaks.kmeans.wss.prop,delta.rand.distpeaks.kmeans.wss.prop,by="row.names"),id.vars = "Row.names",variable.name = "Run",value.name = "Prop"),prob.distpeaks.kmeans.wss,by.x = "Row.names",by.y = "row.names") %>%
      mutate(Cluster=factor(Row.names,levels = paste("k",2:max_k,sep = "")),Run=ifelse(Run=="Obs","obs","rand"),Signif=factor(Signif)) %>%
      select(-Row.names)
    rm(list = c("delta.distpeaks.kmeans.wss.prop","delta.rand.distpeaks.kmeans.wss.prop","prob.distpeaks.kmeans.wss"));gc()
  } else {
    K = max_k
    kopt.peaks.wss.prop <- NULL
    delta.kopt.peaks.wss.prop <- NULL
  }
  
  # Calculation of final kmeans
  cat("Calculating final k-clusters\n")
  peaks.kmeans <- kmeans(peaks_ts.nonzero,centers=K,iter.max = 100,nstart = 1000)
  
  # Adding cluster info to df
  peaks_ts.nonzero.clustered <- merge(data.frame(peaks_ts.nonzero,check.names = F) %>% mutate(peakco=rownames(.)),
                                      data.frame(Cluster=peaks.kmeans$cluster),by.x="peakco",by.y="row.names") %>%
    select(peakco,Cluster,everything()) %>%
    arrange(Cluster)
  
  # Reorders/renumbers the clusters to place more "closing" peak clusters at the top and more "opening" clusters at the bottom
  # These effectively renames all clusters
  clustage.medians <- peaks_ts.nonzero.clustered %>%
    select(-peakco) %>%
    group_by(Cluster) %>% 
    summarize_all(median) %>% 
    select(-Cluster) %>% 
    t() %>% 
    data.frame(.,binage,check.names = F) %>%
    group_by(binage) %>% 
    summarize_all(median) %>% 
    data.frame(.,row.names = "binage",check.names = F) %>%
    t() %>% 
    data.frame(.,row.names = paste("Cluster",rownames(.),sep=""),check.names = F) %>% 
    mutate(mranks=apply(.,1,function(x) mean(diff(rank(x))))) %>%
    mutate(Cluster=row_number()) %>%
    arrange(mranks) %>%
    mutate(cranks=order(mranks))
  
  peaks_ts.nonzero.clustered.sorted <- peaks_ts.nonzero.clustered %>%
    left_join(clustage.medians %>% select(Cluster,cranks),by="Cluster") %>%
    # mutate(porder=apply(peaks_ts.clustered %>% select(-Cluster,-peakco),1,function(x) mean(diff(x)))) %>%
    mutate(porder=apply(peaks_ts.nonzero.clustered %>% select(-Cluster,-peakco),1,function(x) max((x-mean(x))/sd(x)))) %>%
    select(-Cluster) %>%
    rename_(Cluster="cranks") %>%
    arrange(Cluster,-porder) %>%
    select(peakco,Cluster,everything(),-porder)
  
  peaks_ts.nonzero.zscores.clustered.sorted <- peaks_ts.nonzero.clustered.sorted %>%  # z-scores in long format
    select(-Cluster,-peakco) %>%
    apply(.,1,function(x) if (!kmeans_on_zscores) {(x-mean(x))/sd(x)} else {x}) %>%
    t() %>%
    data.frame(Cluster=peaks_ts.nonzero.clustered.sorted$Cluster,.,check.names = F) %>%
    melt(.,id.vars="Cluster",variable.name="aget",value.name="score") %>%
    mutate(aget=as.character(aget),
           age=as.numeric(gsub("[A-z]","",aget)),
           sex=ifelse(grepl("f|m",aget),ifelse(grepl("f",aget),"Females","Males"),"")) %>%
    left_join(data.frame(bp=bincuts) %>% mutate(aget=rownames(.)),by="aget") %>%
    left_join(peaks_ts.nonzero.clustered.sorted %>% group_by(Cluster) %>% summarize(ClusterSize=n()),by="Cluster") %>%
    mutate(Cluster.label=factor(paste("Cluster ",Cluster," (n=",ClusterSize,")",sep=""),levels = unique(paste("Cluster ",Cluster," (n=",ClusterSize,")",sep=""))))
  
  peaks_ts.nonzero.zscores.clustered.sorted.medians <- peaks_ts.nonzero.zscores.clustered.sorted %>%
    group_by(Cluster.label,aget) %>%
    summarize(medians=median(score)) %>%
    dcast(Cluster.label~aget,value.var = "medians") %>%
    data.frame(row.names = "Cluster.label",check.names = F)
  
  peaks_ts.nonzero.zscores.clustered.sorted.medians.hclust <- hclust(dist(peaks_ts.nonzero.zscores.clustered.sorted.medians,method = "maximum"),method = "ward.D2")
  
  # Apply reclustering if using
  if (!is.null(recluster_k)) {
    recluster_k = min(recluster_k,K)
    peaks_ts.nonzero.clustered.sorted <- peaks_ts.nonzero.clustered.sorted %>% 
      left_join(data.frame(Cluster=1:K,reCluster=cutree(peaks_ts.nonzero.zscores.clustered.sorted.medians.hclust,recluster_k)),by="Cluster") %>%
      select(peakco,Cluster,reCluster,everything()) %>%
      arrange(reCluster,Cluster)
    peaks_ts.nonzero.zscores.clustered.sorted <- peaks_ts.nonzero.zscores.clustered.sorted %>% 
      left_join(data.frame(Cluster=1:K,reCluster=cutree(peaks_ts.nonzero.zscores.clustered.sorted.medians.hclust,recluster_k)),by="Cluster") %>%
      left_join(peaks_ts.nonzero.clustered.sorted %>% group_by(reCluster) %>% summarize(reClusterSize=n()),by="reCluster") %>%
      mutate(reCluster.label=factor(paste("reCluster ",reCluster," (n=",reClusterSize,")",sep=""),levels = unique(paste("reCluster ",reCluster," (n=",reClusterSize,")",sep=""))))
    peaks_ts.nonzero.zscores.reclustered.sorted.medians <- peaks_ts.nonzero.zscores.clustered.sorted %>%
      group_by(reCluster.label,aget) %>%
      summarize(medians=median(score)) %>%
      dcast(reCluster.label~aget,value.var = "medians") %>%
      data.frame(row.names = "reCluster.label",check.names = F)
  } else {
    peaks_ts.nonzero.clustered.sorted <- peaks_ts.nonzero.clustered.sorted %>%
      mutate(reCluster=Cluster) %>%
      select(peakco,Cluster,reCluster,everything())
    peaks_ts.nonzero.zscores.clustered.sorted <- peaks_ts.nonzero.zscores.clustered.sorted %>%
      mutate(reCluster=Cluster,
             reClusterSize=ClusterSize,
             reCluster.label=Cluster.label)
    peaks_ts.nonzero.zscores.reclustered.sorted.medians <- peaks_ts.nonzero.zscores.clustered.sorted.medians %>%
      data.frame(.,row.names = sub("^","re",rownames(.)),check.names = F)
  }
  
  # Export of peak-optimal kmean configurations
  saveobj <- c("peaks_ts.nonzero",
               "peaks_ts.nonzero.clustered.sorted",
               "peaks_ts.nonzero.zscores.clustered.sorted",
               "peaks_ts.nonzero.zscores.clustered.sorted.medians",
               "peaks_ts.nonzero.zscores.clustered.sorted.medians.hclust",
               "peaks_ts.nonzero.zscores.reclustered.sorted.medians",
               "peaks.kmeans",
               "K",
               "kopt.peaks.wss.prop",
               "delta.kopt.peaks.wss.prop",
               "nrand",
               "p.rand",
               "recluster_k",
               "max_k",
               "kmeans_on_zscores",
               "saveobj")
  output <- setNames(list(peaks_ts.nonzero,
                          peaks_ts.nonzero.clustered.sorted,
                          peaks_ts.nonzero.zscores.clustered.sorted,
                          peaks_ts.nonzero.zscores.clustered.sorted.medians,
                          peaks_ts.nonzero.zscores.clustered.sorted.medians.hclust,
                          peaks_ts.nonzero.zscores.reclustered.sorted.medians,
                          peaks.kmeans,
                          K,
                          kopt.peaks.wss.prop,
                          delta.kopt.peaks.wss.prop,
                          nrand,
                          p.rand,
                          recluster_k,
                          max_k,
                          kmeans_on_zscores,
                          saveobj),saveobj)
}

###################################################################################################################################################################
# Scripts for KMA applied to ATAC-seq data.
# Three instances are included here, depending on the treatment of sex-specific and sex-agnostic peaks:
# Female (Male) -specific analysis - including all peaks called as trendy in females (males), and peaks that are called as trendy for both sexes
###################################################################################################################################################################

## Load base data and annotations
peakSet="narrow"
load(paste0("./filtered_global.stats_",peakSet,"Peaks_filtered.RData"))
load(paste0("./data/",peakSet,"Peaks/pbmc_aging.summer2017_",peakSet,"Peaks_filtered.RData"))
run.name = "AgexSex"
batch.variable = "batchdate"

###################################################################################################################################################################
## 1 - Female-specific analysis
###################################################################################################################################################################

# Options used in manuscript
nclust=5
K=3

# Female differential peaks
females_fname <- paste0("./data/",peakSet,"Peaks/peak.clusters_selection/",nclust,"clusters/aging.summer2017_",peakSet,"_peaks.kbins_trendy_females_",run.name,"_",batch.variable,".RData")
print("Processing females trendy peaks...")

# Female trendy peaks
kbins_output <- load(paste0("./data/",peakSet,"Peaks/kbins_age_output_females_",run.name,"_",batch.variable,".RData"))
arima_output <- load(paste0("./data/",peakSet,"Peaks/arima_output_females_",run.name,"_",batch.variable,".RData"))
trendypeaks.females <- rownames(predicted.nonzero.females)
binage.females <- setNames(.bincode(evalage,c(0,(bincuts.females)-1,Inf)),evalage)
tsage.females <- sort(unique(age[sex=="F" & age %in% names(binage.females)]))
binage.females.obs <- setNames(.bincode(tsage.females,c(0,(bincuts.females)-1,Inf)),tsage.females)
fitted.nonzero.females <- data.frame(arima.females$fitted.nonzero,check.names = F) %>% select(which(colnames(.) %in% tsage.females))

kbins_atac.fitted.females <- kbins_peaks(peaks_ts=fitted.nonzero.females %>% rename_all(funs(sub("$","f",.))),
                                         binage=setNames(binage.females.obs,sub("$","f",names(binage.females.obs))),
                                         max_k = nclust,
                                         recluster_k = K,
                                         nrand=0,
                                         p.rand=0.05,
                                         seed=1,
                                         kmeans_on_zscores=T,
                                         dataType = "atac",
                                         n.cores=4)
list2env(kbins_atac.fitted.females,envir=Temp)
for (df in Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)]) {
  assign(paste(df,"females",sep="_"),get(df,envir = Temp),pos = -1)
}
save(file = females_fname,list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"females",sep = "_"))

## Export annotated BED files for (fine) clusters (and reclusters)
write.table(atac.allpeaks.anno.pbmc[peaks_ts.nonzero.clustered.sorted_females$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
            file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.trendy_females.bed"),
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")
for (k in 1:recluster_k_females) {
  write.table(atac.allpeaks.anno.pbmc[(peaks_ts.nonzero.clustered.sorted_females %>% filter(reCluster==k))$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
              file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.kbins_trendy_females_cluster",k,".bed"),
              # file = paste0("./data/",peakSet,"Peaks/motif.enrichment/",nclust,"clusters/",peakSet,"Peaks.kbins_trendy_females_cluster",k,".bed"),
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")
}
rm(list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"females",sep = "_"))
rm(list = ls(Temp),envir = Temp)
  
###################################################################################################################################################################
## Male-specific analysis
###################################################################################################################################################################

# Options used in manuscript
nclust=6
K=3

males_fname <- paste0("./data/",peakSet,"Peaks/peak.clusters_selection/",nclust,"clusters/aging.summer2017_",peakSet,"_peaks.kbins_trendy_males_",run.name,"_",batch.variable,".RData")
print("Processing males trendy peaks...")

# Male trendy peaks
kbins_output <- load(paste0("./data/",peakSet,"Peaks/kbins_age_output_males_",run.name,"_",batch.variable,".RData"))
arima_output <- load(paste0("./data/",peakSet,"Peaks/arima_output_males_",run.name,"_",batch.variable,".RData"))
trendypeaks.males <- rownames(predicted.nonzero.males)
binage.males <- setNames(.bincode(evalage,c(0,(bincuts.males)-1,Inf)),evalage)
tsage.males <- sort(unique(age[sex=="M" & age %in% names(binage.males)]))
binage.males.obs <- setNames(.bincode(tsage.males,c(0,(bincuts.males)-1,Inf)),tsage.males)
fitted.nonzero.males <- data.frame(arima.males$fitted.nonzero,check.names = F) %>% select(which(colnames(.) %in% tsage.males))

kbins_atac.fitted.males <- kbins_peaks(peaks_ts=fitted.nonzero.males %>% rename_all(funs(sub("$","m",.))),
                                       binage=setNames(binage.males.obs,sub("$","m",names(binage.males.obs))),
                                       max_k=6,
                                       recluster_k = 3,
                                       nrand=0,
                                       p.rand=0.05,
                                       seed=20090801,
                                       kmeans_on_zscores=T,
                                       dataType = "atac",
                                       n.cores=4)
list2env(kbins_atac.fitted.males,envir=Temp)
for (df in Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)]) {
  assign(paste(df,"males",sep="_"),get(df,envir = Temp),pos = -1)
}
save(file = males_fname,list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"males",sep = "_"))

## Export annotated BED files for (fine) clusters (and reclusters)
write.table(atac.allpeaks.anno.pbmc[peaks_ts.nonzero.clustered.sorted_males$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
            file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.trendy_males.bed"),
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")
for (k in 1:recluster_k_males) {
  write.table(atac.allpeaks.anno.pbmc[(peaks_ts.nonzero.clustered.sorted_males %>% filter(reCluster==k))$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
              file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.kbins_trendy_males_cluster",k,".bed"),
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")
}
rm(list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"males",sep = "_"))
rm(list = ls(Temp),envir = Temp)
  
###################################################################################################################################################################
## F/M-intersection analysis (depends on runs above)
###################################################################################################################################################################

# Options used in manuscripts
nclust=3
K=3

common_fname <- paste0("./data/",peakSet,"Peaks/peak.clusters_selection/",nclust,"clusters/aging.summer2017_",peakSet,"_peaks.kbins_trendy_common_",run.name,"_",batch.variable,".RData")
print("Processing common trendy peaks...")

predicted.common <- data.frame(merge(data.frame(predicted.nonzero.females,check.names = F) %>% rename_all(funs(sub("$","f",.))),
                                     data.frame(predicted.nonzero.males,check.names = F) %>% rename_all(funs(sub("$","m",.))),
                                     by="row.names"),row.names = 1,check.names = F)
trendypeaks.common <- rownames(predicted.common)
fitted.common <- data.frame(merge(data.frame(fitted.nonzero.females,check.names = F) %>% rename_all(funs(sub("$","f",.))),
                                  data.frame(fitted.nonzero.males,check.names = F) %>% rename_all(funs(sub("$","m",.))),
                                  by="row.names"),row.names = 1,check.names = F)

kbins_atac.fitted.common <- kbins_peaks(peaks_ts=fitted.common,
                                        binage=c(setNames(binage.females.obs,sub("$","f",names(binage.females.obs))),setNames(binage.males.obs,sub("$","m",names(binage.males.obs)))),
                                        max_k=nclust,
                                        recluster_k = K,
                                        nrand=0,
                                        p.rand=0.05,
                                        seed=20090801,
                                        kmeans_on_zscores=T,
                                        dataType = "atac",
                                        n.cores=4)
list2env(kbins_atac.fitted.common,envir=Temp)
for (df in Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)]) {
  assign(paste(df,"common",sep="_"),get(df,envir = Temp),pos = -1)
}
save(file = common_fname,list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"common",sep = "_"))

## Export annotated BED files for (fine) clusters (and reclusters)
write.table(atac.allpeaks.anno.pbmc[peaks_ts.nonzero.clustered.sorted_common$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
            file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.trendy_common.bed"),
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")
for (k in 1:recluster_k_common) {
  write.table(atac.allpeaks.anno.pbmc[(peaks_ts.nonzero.clustered.sorted_common %>% filter(reCluster==k))$peakco,c("chr","start","end","GeneName","RE_specificity_fine")],
              file = paste0("./data/",peakSet,"Peaks/motif.enrichment/final/",peakSet,"Peaks.kbins_trendy_common_cluster",k,".bed"),
              quote = F,
              row.names = F,
              col.names = F,
              sep = "\t")
}
rm(list = paste(Temp$saveobj[grep("saveobj",Temp$saveobj,invert = T)],"common",sep = "_"))
rm(list = ls(Temp),envir = Temp)
