# Reads motifs (MEME format) and motif comparison output (from TOMTOM), finds clusters of motifs and outputs finalized RData file

library(igraph)
library(tibble)
library(dplyr)
# Results in separate environment

motif.env = new.env()
motif.env$knownmotif.path = "./data-repo/knownMotifs_meme.txt"
known.motifs_aliases <- read.table("./data-repo/cluster_aliases.txt",header = T,quote = '"',sep = "\t",stringsAsFactors = F) # Pre-annotated clusters (matches clusters/components calculated below)
known.motifs <- with(motif.env,readLines(con = knownmotif.path,n = -1))
known.motifs_idx <- grep("^letter\\-probability",known.motifs)
known.motifs_len <- as.numeric(sapply(known.motifs_idx, function(idx) sub(" .*","",sub(".*w\\= ","",known.motifs[idx]))))
known.motifs_names <- sapply(known.motifs_idx, function(idx) sub(".* ","",known.motifs[idx-2]))
known.motifs_consensus <- sapply(known.motifs_idx, function(idx) sub(" .*","",sub("MOTIF ","",known.motifs[idx-2])))
motif.env$knownMotifsIDs <- data.frame(MotifName=known.motifs_names,Consensus=known.motifs_consensus,stringsAsFactors = F)

# Extract and format known motifs
motif.env$knownMotifs <- lapply(setNames(seq(length(known.motifs_idx)),known.motifs_names), function(id) {
  idx = known.motifs_idx[id]
  len = known.motifs_len[id]
  apply(matrix(sapply(strsplit(known.motifs[(idx+1):(idx+len)],split = "\t"),as.numeric),nrow = 4,byrow = F,dimnames = list(c("A","C","G","T"))),2,function(x) x/sum(x))
})

# Collect family information, tagging redundant TFs
Evalue.thresh = 5
qvalue.thresh = 0.00001
motif.env$knownMotifsComparison <- read.table("./data-repo/tomtom_out/tomtom.tsv",header = T,sep = "\t",quote = "",stringsAsFactors = F)

known.motifs.network <- do.call(rbind,lapply(known.motifs_consensus,function(motif) {
  motif.env$knownMotifsComparison %>% 
    filter(Query_ID==motif) %>%
    filter(q.value<=qvalue.thresh,E.value<=Evalue.thresh) %>%
    select(Query_ID,Target_ID)
}))
known.motifs.clusters <- components(as.undirected(simplify(graph.data.frame(known.motifs.network))),"strong")

motif.env$knownMotifsClusterMembership <- data.frame(Cluster=known.motifs.clusters$membership) %>%
  rownames_to_column("Motif") %>%
  left_join(known.motifs_aliases,by="Cluster") %>%
  arrange(Cluster,Motif,Alias)
motif.env$knownMotifsClusters <- motif.env$knownMotifsClusterMembership %>%
  distinct(Cluster,Alias,N,TFNames)

known.motifs.detail <- data.frame(Cluster=known.motifs.clusters$membership) %>%   # to match clusters/motifs to actual TF names
  rownames_to_column("Consensus") %>%
  inner_join(motif.env$knownMotifsIDs,by="Consensus")

rm(list = ls(pattern = "known.motifs*"))
save(motif.env,file = "./data-repo/knownMotifsData.RData")
detach("package:igraph", unload=TRUE)