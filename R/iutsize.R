iutsize <- function(p, n, alpha=0.1, sim=1e6){
  
  i <- count <- 0
  
  q <- qf(1 - alpha, df1=p, df2=n - p) / ((n - p) / (p * (n - 1)))
  
  while(i < sim){
    
    if(rf(1, df1=1, df2=n - 1) > q){
      count <- count + 1
    }
    
    i <- i + 1
    
  }
  
  sum(count) / (2 * sim)
  
}