confints <- function(dat, method, alpha=0.1, steps=100){
  
  method <- match.arg(method, choices=c("expanded", "fixseq", "tost"))
  
  if(method=="expanded"){
    
    n <- nrow(dat)
    
    est <- colMeans(dat)
    sd <- apply(dat, 2, sd)
    
    ci <- matrix(c(est - sd * qt(1 - alpha, n - 1) / sqrt(n),
                   est + sd * qt(1 - alpha, n - 1) / sqrt(n)), ncol(dat))
    ci[, 1] <- ifelse(ci[, 1] > 0, 0, ci[, 1])
    ci[, 2] <- ifelse(ci[, 2] < 0, 0, ci[, 2])
    colnames(ci) <- c("lower", "upper")
    
  }
  
  if(method=="fixseq"){
    
    stop("Not yet.")
    
  }
  
  if(method=="tost"){
    
    n <- nrow(dat)
    
    est <- colMeans(dat)
    sd <- apply(dat, 2, sd)
    
    ci <- matrix(c(est - sd * qt(1 - alpha, n - 1) / sqrt(n),
                   est + sd * qt(1 - alpha, n - 1) / sqrt(n)), ncol(dat))
    colnames(ci) <- c("lower", "upper")
    
  }
  
  Out <- cbind(est, ci)
  rownames(Out) <- colnames(dat)
  colnames(Out) <- c("estimate", "lower", "upper")
  
  return(Out)
  
}