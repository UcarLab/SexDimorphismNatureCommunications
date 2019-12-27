# Main script used to fit time series ARIMA models to (batch-adjusted) ATAC-seq data.
# Adjusted data is derived as a by-product of differential analyses.

source("/path/to/arima.selection.R")

ts_atac <- function(adj.data=NULL,
                    metadataFactors=NULL, # sex, age, sampid, etc
                    evalage=NULL,  # Range of ages where time series models are to be evaluated for plotting, clustering, etc. For example, all ages between the min and max age, or all ages included in the sample
                    arima.use.boxcox=FALSE,
                    arima.use.portmanteau=TRUE,
                    arima.loess.span=0.75,
                    randtest.max_k=2,
                    randtest.n=1000,
                    randtest.seed=1,
                    randtest.use.zscores=FALSE,
                    output.tmp="./",
                    dataType=c("atac","rna"),
                    n.cores=1) {
  
  cat("Starting run:\n")
  
  cat("Fitting ARIMA model to peaks: finding trendy peaks, computing fitted and predicted data\n")
  predicted.nonzero = list()
  if (any("F" %in% metadataFactors[colnames(adj.data),]$sex)) {
    cat("ARIMA fitting for females...\n============================\n")
    arima.females <- arima.selection(adj.data=adj.data,
                                     metadata=metadataFactors,
                                     selected.sex="F",
                                     evalage=evalage,
                                     use.boxcox=arima.use.boxcox,
                                     use.portmanteau=arima.use.portmanteau,
                                     loess.span=arima.loess.span,
                                     dataType=dataType[[1]],
                                     n.cores=n.cores)
    save(arima.females,file = paste(output.tmp,"/","tmp_arima.females_",globalPars$run.name,".RData",sep=""))
    predicted.females.nonzero <- arima.females$predicted.nonzero
    rm("arima.females")
    gc()
    predicted.nonzero$females <- predicted.females.nonzero
  }
  if (any("M" %in% metadataFactors[colnames(adj.data),]$sex)) {
    cat("ARIMA fitting for males...\n============================\n")
    arima.males <- arima.selection(adj.data=adj.data,
                                   metadata=metadataFactors,
                                   selected.sex="M",
                                   evalage=evalage,
                                   use.boxcox=arima.use.boxcox,
                                   use.portmanteau=arima.use.portmanteau,
                                   loess.span=arima.loess.span,
                                   dataType=dataType[[1]],
                                   n.cores=n.cores)
    save(arima.males,file = paste(output.tmp,"/","tmp_arima.males_",globalPars$run.name,".RData",sep=""))
    predicted.males.nonzero <- arima.males$predicted.nonzero
    rm("arima.males")
    gc()
    predicted.nonzero$males <- predicted.males.nonzero
  }
  
  if (any("F" %in% metadataFactors[colnames(adj.data),]$sex)) {
    load(paste(output.tmp,"/","tmp_arima.females_",globalPars$run.name,".RData",sep=""))
    load(paste(output.tmp,"/","tmp_ts.age.kbins.females_",globalPars$run.name,".RData",sep=""))
  }
  if (any("M" %in% metadataFactors[colnames(adj.data),]$sex)) {
    load(paste(output.tmp,"/","tmp_arima.males_",globalPars$run.name,".RData",sep=""))
    load(paste(output.tmp,"/","tmp_ts.age.kbins.males_",globalPars$run.name,".RData",sep=""))
  }
  
  output = list()
  if (any("F" %in% metadataFactors[colnames(adj.data),]$sex)) {
    output$arima.females=arima.females
    output$ts.age.kbins.females=ts.age.kbins.females
    output$predicted.nonzero.females=arima.females$predicted.nonzero
  }
  if (any("M" %in% metadataFactors[colnames(adj.data),]$sex)) {
    output$arima.males=arima.males
    output$ts.age.kbins.males=ts.age.kbins.males
    output$predicted.nonzero.males=arima.males$predicted.nonzero
  }
  
  return(output)
}