
# Uses extensive auto-arima to select peaks/transcripts showing evidence of non-negligible trend based on adjusted peak/transcript data,
# then uses LOESS regression to fit a continuous time series to each peak

arima.selection <- function(adj.data=NULL,
                            metadata=NULL,
                            selected.sex=NULL,
                            evalage=NULL,
                            use.boxcox=FALSE,
                            use.portmanteau=TRUE,
                            loess.span=0.75,
                            dataType=c("atac","rna"),
                            n.cores=21) {
  
  require(zoo)
  require(forecast)
  require(reshape2)
  require(parallel)
  require(dplyr)
  
  age=with(metadata,setNames(age,sampid))
  sex=with(metadata,setNames(sex,sampid))

  cat("Building time series object\n")
  age.sexsel <- age[names(sex[sex==selected.sex])]
  age.sexsel.sorted <- age.sexsel[order(age.sexsel)]
  adj.data.sexsel <- adj.data[,names(sex[sex==selected.sex])]
  
  selpeaks.sexsel <- t(adj.data.sexsel)
  age.counts <- as.data.frame(table(age.sexsel),stringsAsFactors = F)
  selpeaks.sexsel.unique <- do.call(rbind,lapply(setNames(as.list(seq_len(nrow(age.counts))),age.counts$age.sexsel),function(n) {
    if (age.counts[n,"Freq"]==1) {
      estage <- selpeaks.sexsel[age.sexsel==age.counts[n,"age.sexsel"],]
    } else {
      estage <- colMeans(selpeaks.sexsel[age.sexsel==age.counts[n,"age.sexsel"],])
    }
    return(estage)}))
  age.sexsel.unique <- as.numeric(age.counts$age.sexsel)

  # Create zoo/ts objects, average data for same ages, interpolate unsampled ages
  zoo.sexsel <- zoo(selpeaks.sexsel.unique,order.by=age.sexsel.unique,frequency = 1)
  ts.sexsel <- mclapply(zoo.sexsel,ts,mc.cores = n.cores)
  # plot.ts(ts.sexsel[[1]],type="b") # Visualization of a single TS
  if (use.boxcox) {
    cat("Calculating Box-Cox Lambda using the log-likelihood method\n")
    bclambda.sexsel <- sapply(ts.sexsel,BoxCox.lambda,method="loglik",upper=10)
  } else {
    bclambda.sexsel = NULL
  }
  gc()
  
  cat("Fitting auto-arima models on",n.cores,"CPU cores\n")
  # Builds an ARIMA model to test for trends in each time series, then selects those for which the trend is non-negligible (stationarity cannot be rejected AND no AR or MA terms are chosen)
  arima.sexsel <- mclapply(setNames(seq_len(length(ts.sexsel)),names(ts.sexsel)),function(n) auto.arima(ts.sexsel[[n]],seasonal = F,allowmean = T,allowdrift = T,lambda = bclambda.sexsel), mc.cores = n.cores)
  gc()
  if (use.portmanteau) {
    # portmanteau tests on RESIDUALS (Ljung-Box) - models with autocorrelated residuals are removed as well
    cat("Computing portmanteau tests on ARIMA model residuals using the Ljung-Box algorithm, on",n.cores,"CPU cores\n")
    portmanteau.sexsel <- unlist(mclapply(arima.sexsel,function(x) {
      df=sum(x$arma[1:2])
      pmin <- Box.test(x$residuals,lag=20,type="Ljung-Box",fitdf = df)$p.value
      return(pmin)
    },mc.cores = n.cores))
    cat("Proceeding to filter out",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0)) & portmanteau.sexsel>0.05
  } else {
    cat("Proceeding to filter out",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0))
    portmanteau.sexsel = NULL
  }
  gc()

  arima.sexsel.nonzero <- arima.sexsel[arima.sexsel.select]
  ts.sexsel.nonzero <- ts.sexsel[arima.sexsel.select]
  fitted.sexsel.nonzero <- t(as.data.frame(sapply(arima.sexsel.nonzero,`[[`,"fitted"),row.names = as.integer(attr(arima.sexsel.nonzero[[1]]$x,"index"))))
  cat("Out of",length(arima.sexsel),ifelse(dataType[[1]]=="atac","peaks,","transcripts,"),length(ts.sexsel.nonzero),"trendy and well-fitted",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"were selected for further analysis\n")
  gc()
  
  ## Uses local fit (LOESS) on sex-union nonzero trends to predict chromatin values at each age in the span, for sex comparisons
  cat("Computing predicting",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"using LOESS interpolation over ARIMA-fitted data on",n.cores,"cores\n")
  loess.sexsel.nonzero <- mclapply(arima.sexsel.nonzero,function(A) loess(as.vector(A$fitted)~attr(A$x,"index"),span = loess.span,family = "symmetric"),mc.cores = n.cores)
  gc()
  predicted.sexsel.nonzero <- t(as.data.frame(mclapply(loess.sexsel.nonzero,predict,evalage,mc.cores = n.cores),row.names = evalage))
  gc()
  
  output <- list(ts.complete=ts.sexsel,
                 arima.complete=arima.sexsel,
                 arima.select=arima.sexsel.select,
                 loess.nonzero=loess.sexsel.nonzero,
                 fitted.nonzero=fitted.sexsel.nonzero,
                 predicted.nonzero=predicted.sexsel.nonzero,
                 BoxCox_lambda=bclambda.sexsel,
                 portmanteau.pvalues=portmanteau.sexsel)
  
  cat("Done.\n")
  return(output)
}