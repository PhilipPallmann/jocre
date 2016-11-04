translate <- function(level=0.95, ddf, direction){
  
  dir <- match.arg(direction, choices=c("ci2cr", "cr2ci"))
  
  if(dir=="ci2cr"){
    
    if(ddf==0){
      out <- pchisq(qchisq(level, df=1), df=2)
    }else{
      out <- pf(0.5 * qf(level, df1=1, df2=ddf), df1=2, df2=ddf)
    }
    
  }else{
    
    if(ddf==0){
      out <- pchisq(qchisq(level, df=2), df=1)
    }else{
      out <- pf(2 * qf(level, df1=2, df2=ddf), df1=1, df2=ddf)
    }
  }
  
  return(out)
  
}