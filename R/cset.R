cset <- function(dat, method, alpha=0.1, steps=100, TsengBrownA=1, TsengBrownB=1){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  est <- matrix(colMeans(dat), p)
  poolvar <- var(as.vector(as.matrix(dat)))
  cov <- cov(dat)
  solved <- solve(cov)
  #s2 <- poolvar/n
  #s <- sqrt(poolvar/n)
  
  method <- match.arg(method, choices=c("bootkern", "emp.bayes", "hotelling", "limacon.asy", "limacon.fin",
                                        "standard.cor", "standard.ind", "tseng", "tseng.brown"))
  
  if(method=="bootkern"){
    
    stop("Not implemented (yet).")
    
  }
  
  if(method=="emp.bayes"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      a <- (p - 2) * (df / (df + 2))
      JSfactor <- (1 - (a * (poolvar/n)) / sum(est^2)) # or just poolvar???
      
      if(JSfactor < 0){
        JSplus <- est
      }else{
        JSplus <- JSfactor * est
      }
      
      cond <- (sum(est^2) / (poolvar/n)) < (p * qf(p=1 - alpha, df1=p, df2=df)) # or just poolvar???
      
      if(cond==TRUE){
        vE2 <- (1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))) *
          (p * qf(p=1 - alpha, df1=p, df2=df) - p * log(1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))))
      }else{
        vE2 <- (1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))) *
          (p * qf(p=1 - alpha, df1=p, df2=df) - p * log(1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))))
      }
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        sqrt(sum((theta - JSplus)^2)) < (sqrt(poolvar/n) * sqrt(vE2))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="hotelling"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      rhs <- qf(p=1 - alpha, df1=p, df2=df - p + 1) * p * df / (df - p + 1)
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        (n * t(est - theta) %*% solved %*% (est - theta)) < rhs 
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="limacon.asy"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% solved %*% est) / sqrt(t(theta) %*% solved %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solved %*% theta) * sqrt(n)) 
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="limacon.fin"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="standard.cor"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        t(est - theta) %*% solved %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=df))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="standard.ind"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < (poolvar/n * p * qf(p=1 - alpha, df1=p, df2=df))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="tseng"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      lhs <- sqrt(sum(est^2))^2 / (p * poolvar/n)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        lhs > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / (poolvar/n)))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
      if(nrow(crFinal)==0){
        ciFinal <- matrix(rep(Inf, 2 * p), 2)
        break
      }
      
      if(ciFinal[1, 1]==Inf | ciFinal[1, 1]==-Inf){
        ciFinal <- matrix(rep(0, 2 * p), 2)
        break
      }
      
    }
    
  }
  
  if(method=="tseng.brown"){
    
    ciFinal <- cbind(est - 2 * poolvar, est + 2 * poolvar)
    
    searchwidth <- 1
    
    while(min(abs(ciFinal[, 1] - (est - 2 * poolvar))) < 0.001 | min(abs(ciFinal[, 2] - (est + 2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
          qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
      })
      
      crFinal <- cbind(grid, findcr)[findcr==1, ]
      ciFinal <- t(apply(crFinal[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
            
      if(nrow(crFinal)==0){
        ciFinal <- matrix(rep(Inf, 2 * p), 2)
        break
      }
      
      if(ciFinal[1, 1]==Inf | ciFinal[1, 1]==-Inf){
        ciFinal <- matrix(rep(0, 2 * p), 2)
        break
      }
      
    }
    
  }
  
  Out <- list()
  
  Out$cr <- crFinal
  Out$ci <- ciFinal
  Out$n <- n
  Out$p <- p
  Out$est <- est
  Out$cov <- cov
  Out$poolvar <- poolvar
  Out$df <- df
  Out$dat <- dat
  Out$method <- method
  Out$alpha <- alpha
  Out$steps <- steps
  
  class(Out) <- "JOC"
  
  return(Out)
  
}
