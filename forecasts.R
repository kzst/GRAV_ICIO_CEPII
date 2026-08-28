forecasts <-function(fit,h=10,level=95,bayesian=FALSE){
  forecasts<-list()
  for (i in 1:length(fit)){
    if (bayesian){
      f<-bayesforecast::forecast(fit[[i]],h=h,probs=level/100)
      #f<-as.data.frame(f)
      #colnames(f)
      #colnames(f)<-c("Point Forecast",paste("Lo",level),paste("Hi",level))
      forecasts[[i]]<-f
    }else{
      forecasts[[i]]<-forecast::forecast(fit[[i]],h=h,level=level)
    }
    
  }
  return(forecasts)
}
astibble <-function(fc,index=2021:2030){
  astibble<-list()
  for (i in 1:length(fc)){
    astibble[[i]] <- fc[[i]] %>% 
      as_tibble() %>% 
      mutate(index = index)
  }
  return(astibble)
}