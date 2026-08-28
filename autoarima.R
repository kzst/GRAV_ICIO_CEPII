autoarima <-function(data,start=1995,end=2020,bayesian=FALSE){
  fit<-list()
  for (i in 1:ncol(data)){
    TS<-ts(data[,i],start=start,end=end)
    if (bayesian){
      fit[[i]]<-bayesforecast::auto.sarima(TS,seasonal = FALSE)      
    }else{
      fit[[i]]<-forecast::auto.arima(TS)
    }
  }
  return(fit)
}