confset <- function(dat, method, alpha=0.1, steps=100, TsengBrownA=1, TsengBrownB=1){
  
  method <- match.arg(method, choices=c("bootkern", "emp.bayes", "hotelling", "limacon.asy", "limacon.fin",
                                        "standard.cor", "standard.ind", "tseng", "tseng.brown")) 
  
  if(method=="bootkern"){
    
    stop("Not yet.")
    
  }
  
  if(method=="emp.bayes"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    a <- (p - 2) * (df / (df + 2))
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    s <- sqrt(poolvar / n)
    
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
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrCas <- apply(grid, 1, function(x){
      theta <- matrix(x, p)
      sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
    })
    
    crCas <- cbind(grid, findcrCas)[findcrCas==1, ]
    
    Cas0 <- t(apply(crCas[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Cas0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Cas0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrCas2 <- apply(grid2, 1, function(x){
        theta <- matrix(x, p)
        sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
      })
      
      crCas2 <- cbind(grid2, findcrCas2)[findcrCas2==1, ]
      
      if(min(crCas2[, 1])==min(grid2[, 1]) | max(crCas2[, 1])==max(grid2[, 1]) |
           min(crCas2[, 2])==min(grid2[, 2]) | max(crCas2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrCas3 <- apply(grid3, 1, function(x){
          theta <- matrix(x, p)
          sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
        })
        
        crCas3 <- cbind(grid3, findcrCas3)[findcrCas3==1, ]
        
        crCas2 <- crCas3
        
      }
      
      Cas0 <- t(apply(crCas2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Cas0[1, 1] - stepwidth, Cas0[1, 1], length.out=steps),
                           seq(Cas0[2, 1] - stepwidth, Cas0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Cas0[1, 2], Cas0[1, 2] + stepwidth, length.out=steps),
                           seq(Cas0[2, 1], Cas0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Cas0[1, 1] - stepwidth, Cas0[1, 2], length.out=steps),
                           seq(Cas0[2, 1] - stepwidth, Cas0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Cas0[1, 1], Cas0[1, 2] + stepwidth, length.out=steps),
                           seq(Cas0[2, 2], Cas0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrCas <- apply(Grid, 1, function(x){
        theta <- matrix(x, 2)
        sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
      })
      
      CrCas <- cbind(Grid, FindcrCas)[FindcrCas==1, ]
      
      Cas <- t(apply(CrCas[, -(p + 1)], 2, range))
      
    }else{
      
      Cas <- Cas0
      
    }
    
    CasOut <- cbind(JSplus, Cas)
    rownames(CasOut) <- colnames(dat)
    colnames(CasOut) <- c("estimate", "lower", "upper")
    
    return(CasOut)
    
  }
  
  if(method=="hotelling"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrHot <- apply(grid, 1, function(x){
      theta <- matrix(x, p)
      (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
        qf(p=1 - alpha, df1=p, df2=df - p + 1) * p * df / (df - p + 1)
    })
    
    crHot <- cbind(grid, findcrHot)[findcrHot==1, ]
    
    Hot0 <- t(apply(crHot[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Hot0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Hot0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrHot2 <- apply(grid2, 1, function(x){
        theta <- matrix(x, p)
        (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
          qf(p=1 - alpha, df1=p, df2=df - p + 1) * p * df / (df - p + 1)
      })
      
      crHot2 <- cbind(grid2, findcrHot2)[findcrHot2==1, ]
      
      if(min(crHot2[, 1])==min(grid2[, 1]) | max(crHot2[, 1])==max(grid2[, 1]) |
           min(crHot2[, 2])==min(grid2[, 2]) | max(crHot2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrHot3 <- apply(grid3, 1, function(x){
          theta <- matrix(x, p)
          (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
            qf(p=1 - alpha, df1=p, df2=df - p + 1) * p * df / (df - p + 1)
        })
        
        crHot3 <- cbind(grid3, findcrHot3)[findcrHot3==1, ]
        
        crHot2 <- crHot3
        
      }
      
      Hot0 <- t(apply(crHot2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Hot0[1, 1] - stepwidth, Hot0[1, 1], length.out=steps),
                           seq(Hot0[2, 1] - stepwidth, Hot0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Hot0[1, 2], Hot0[1, 2] + stepwidth, length.out=steps),
                           seq(Hot0[2, 1], Hot0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Hot0[1, 1] - stepwidth, Hot0[1, 2], length.out=steps),
                           seq(Hot0[2, 1] - stepwidth, Hot0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Hot0[1, 1], Hot0[1, 2] + stepwidth, length.out=steps),
                           seq(Hot0[2, 2], Hot0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrHot <- apply(Grid, 1, function(x){
        theta <- matrix(x, 2)
        (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
          qf(p=1 - alpha, df1=2, df2=df - p + 1) * p * df / (df - p + 1)
      })
      
      CrHot <- cbind(Grid, FindcrHot)[FindcrHot==1, ]
      
      Hot <- t(apply(CrHot[, -(p + 1)], 2, range))
      
    }else{
      
      Hot <- Hot0
      
    }
    
    HotOut <- cbind(est, Hot)
    rownames(HotOut) <- colnames(dat)
    colnames(HotOut) <- c("estimate", "lower", "upper")
    
    return(HotOut)
    
  }
  
  if(method=="limacon.asy"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrLim <- apply(grid, 1, function(x){
      theta <- matrix(x, p)
      ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
        (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
    })
    
    crLim <- cbind(grid, findcrLim)[findcrLim==1, ]
    
    Lim0 <- t(apply(crLim[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Lim0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Lim0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrLim2 <- apply(grid2, 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
      })
      
      crLim2 <- cbind(grid2, findcrLim2)[findcrLim2==1, ]
      
      if(min(crLim2[, 1])==min(grid2[, 1]) | max(crLim2[, 1])==max(grid2[, 1]) |
           min(crLim2[, 2])==min(grid2[, 2]) | max(crLim2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrLim3 <- apply(grid3, 1, function(x){
          theta <- matrix(x, p)
          ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
            (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
        })
        
        crLim3 <- cbind(grid3, findcrLim3)[findcrLim3==1, ]
        
        crLim2 <- crLim3
        
      }
      
      Lim0 <- t(apply(crLim2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 1], length.out=steps),
                           seq(Lim0[2, 1] - stepwidth, Lim0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Lim0[1, 2], Lim0[1, 2] + stepwidth, length.out=steps),
                           seq(Lim0[2, 1], Lim0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 2], length.out=steps),
                           seq(Lim0[2, 1] - stepwidth, Lim0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Lim0[1, 1], Lim0[1, 2] + stepwidth, length.out=steps),
                           seq(Lim0[2, 2], Lim0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrLim <- apply(Grid, 1, function(x){
        theta <- matrix(x, 2)
        ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
      })
      
      CrLim <- cbind(Grid, FindcrLim)[FindcrLim==1, ]
      
      Lim <- t(apply(CrLim[, -(p + 1)], 2, range))
      
    }else{
      
      Lim <- Lim0
      
    }
    
    LimOut <- cbind(est, Lim)
    rownames(LimOut) <- colnames(dat)
    colnames(LimOut) <- c("estimate", "lower", "upper")
    
    return(LimOut)
    
  }
  
  if(method=="limacon.fin"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrLim <- apply(grid, 1, function(x){
      theta <- matrix(x, p)
      ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
        ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
    })
    
    crLim <- cbind(grid, findcrLim)[findcrLim==1, ]
    
    Lim0 <- t(apply(crLim[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Lim0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Lim0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrLim2 <- apply(grid2, 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      crLim2 <- cbind(grid2, findcrLim2)[findcrLim2==1, ]
      
      if(min(crLim2[, 1])==min(grid2[, 1]) | max(crLim2[, 1])==max(grid2[, 1]) |
           min(crLim2[, 2])==min(grid2[, 2]) | max(crLim2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrLim3 <- apply(grid3, 1, function(x){
          theta <- matrix(x, p)
          ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
            ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
        })
        
        crLim3 <- cbind(grid3, findcrLim3)[findcrLim3==1, ]
        
        crLim2 <- crLim3
        
      }
      
      Lim0 <- t(apply(crLim2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 1], length.out=steps),
                           seq(Lim0[2, 1] - stepwidth, Lim0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Lim0[1, 2], Lim0[1, 2] + stepwidth, length.out=steps),
                           seq(Lim0[2, 1], Lim0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 2], length.out=steps),
                           seq(Lim0[2, 1] - stepwidth, Lim0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Lim0[1, 1], Lim0[1, 2] + stepwidth, length.out=steps),
                           seq(Lim0[2, 2], Lim0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrLim <- apply(Grid, 1, function(x){
        theta <- matrix(x, 2)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      CrLim <- cbind(Grid, FindcrLim)[FindcrLim==1, ]
      
      Lim <- t(apply(CrLim[, -(p + 1)], 2, range))
      
    }else{
      
      Lim <- Lim0
      
    }
    
    LimOut <- cbind(est, Lim)
    rownames(LimOut) <- colnames(dat)
    colnames(LimOut) <- c("estimate", "lower", "upper")
    
    return(LimOut)
    
  }
  
  if(method=="standard.cor"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    LambdaInv <- solve(cov(dat))
    poolvar <- var(as.vector(as.matrix(dat)))
    #s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrUsu <- apply(grid, 1, function(x){
      theta <- x
      #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
      t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
    })
    
    crUsu <- cbind(grid, findcrUsu)[findcrUsu==1, ]
    
    Usu0 <- t(apply(crUsu[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Usu0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Usu0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrUsu2 <- apply(grid2, 1, function(x){
        theta <- x
        #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
        t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
      })
      
      crUsu2 <- cbind(grid2, findcrUsu2)[findcrUsu2==1, ]
      
      if(min(crUsu2[, 1])==min(grid2[, 1]) | max(crUsu2[, 1])==max(grid2[, 1]) |
           min(crUsu2[, 2])==min(grid2[, 2]) | max(crUsu2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrUsu3 <- apply(grid3, 1, function(x){
          theta <- x
          #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
          t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
        })
        
        crUsu3 <- cbind(grid3, findcrUsu3)[findcrUsu3==1, ]
        
        crUsu2 <- crUsu3
        
      }
      
      Usu0 <- t(apply(crUsu2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 1], length.out=steps),
                           seq(Usu0[2, 1] - stepwidth, Usu0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Usu0[1, 2], Usu0[1, 2] + stepwidth, length.out=steps),
                           seq(Usu0[2, 1], Usu0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 2], length.out=steps),
                           seq(Usu0[2, 1] - stepwidth, Usu0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Usu0[1, 1], Usu0[1, 2] + stepwidth, length.out=steps),
                           seq(Usu0[2, 2], Usu0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrUsu <- apply(Grid, 1, function(x){
        theta <- x
        #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
        t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
      })
      
      CrUsu <- cbind(Grid, FindcrUsu)[FindcrUsu==1, ]
      
      Usu <- t(apply(CrUsu[, -(p + 1)], 2, range))
      
    }else{
      
      Usu <- Usu0
      
    }
    
    UsuOut <- cbind(est, Usu)
    rownames(UsuOut) <- colnames(dat)
    colnames(UsuOut) <- c("estimate", "lower", "upper")
    
    return(UsuOut)
    
  }
  
  if(method=="standard.ind"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrUsu <- apply(grid, 1, function(x){
      theta <- x
      sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
    })
    
    crUsu <- cbind(grid, findcrUsu)[findcrUsu==1, ]
    
    Usu0 <- t(apply(crUsu[, -(p + 1)], 2, range, na.rm=TRUE))
    
    if(min(abs(Usu0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Usu0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrUsu2 <- apply(grid2, 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
      })
      
      crUsu2 <- cbind(grid2, findcrUsu2)[findcrUsu2==1, ]
      
      if(min(crUsu2[, 1])==min(grid2[, 1]) | max(crUsu2[, 1])==max(grid2[, 1]) |
           min(crUsu2[, 2])==min(grid2[, 2]) | max(crUsu2[, 2])==max(grid2[, 2])){
        
        togrid3 <- list()
        
        for(i in 1:p){
          togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
        }
        
        grid3 <- expand.grid(togrid3)
        
        findcrUsu3 <- apply(grid3, 1, function(x){
          theta <- x
          sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
        })
        
        crUsu3 <- cbind(grid3, findcrUsu3)[findcrUsu3==1, ]
        
        crUsu2 <- crUsu3
        
      }
      
      Usu0 <- t(apply(crUsu2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 1], length.out=steps),
                           seq(Usu0[2, 1] - stepwidth, Usu0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Usu0[1, 2], Usu0[1, 2] + stepwidth, length.out=steps),
                           seq(Usu0[2, 1], Usu0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 2], length.out=steps),
                           seq(Usu0[2, 1] - stepwidth, Usu0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Usu0[1, 1], Usu0[1, 2] + stepwidth, length.out=steps),
                           seq(Usu0[2, 2], Usu0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrUsu <- apply(Grid, 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
      })
      
      CrUsu <- cbind(Grid, FindcrUsu)[FindcrUsu==1, ]
      
      Usu <- t(apply(CrUsu[, -(p + 1)], 2, range))
      
    }else{
      
      Usu <- Usu0
      
    }
    
    UsuOut <- cbind(est, Usu)
    rownames(UsuOut) <- colnames(dat)
    colnames(UsuOut) <- c("estimate", "lower", "upper")
    
    return(UsuOut)
    
  }
  
  if(method=="tseng"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    part1 <- sqrt(sum(est^2))^2 / (p * s2)
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrTse <- apply(grid, 1, function(x){
      theta <- x
      part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
    })
    
    crTse <- cbind(grid, findcrTse)[findcrTse==1, ]
    
    if(nrow(crTse)==0){
      Tse0 <- matrix(rep(Inf, 2 * p), 2)
    }else{
      Tse0 <- t(apply(crTse[, -(p + 1)], 2, range))
    }
    
    if(Tse0[1, 1]==Inf | Tse0[1, 1]==-Inf){
      
      Tse <- matrix(rep(0, 2 * p), 2)
      
    }else{
      
      if(min(abs(Tse0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Tse0[, 2] - est - 2 * poolvar)) < 0.001){
        
        togrid2 <- list()
        
        for(i in 1:p){
          togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
        }
        
        grid2 <- expand.grid(togrid2)
        
        findcrTse2 <- apply(grid2, 1, function(x){
          theta <- x
          part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
        })
        
        crTse2 <- cbind(grid2, findcrTse2)[findcrTse2==1, ]
        
        if(min(crTse2[, 1])==min(grid2[, 1]) | max(crTse2[, 1])==max(grid2[, 1]) |
             min(crTse2[, 2])==min(grid2[, 2]) | max(crTse2[, 2])==max(grid2[, 2])){
          
          togrid3 <- list()
          
          for(i in 1:p){
            togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
          }
          
          grid3 <- expand.grid(togrid3)
          
          findcrTse3 <- apply(grid3, 1, function(x){
            theta <- x
            part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
          })
          
          crTse3 <- cbind(grid3, findcrTse3)[findcrTse3==1, ]
          
          crTse2 <- crTse3
          
        }
        
        Tse0 <- t(apply(crTse2[, -(p + 1)], 2, range))
        
      }
      
      stepwidth <- 4 * poolvar / steps
      
      if(p==2){
        
        GridA <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 1], length.out=steps),
                             seq(Tse0[2, 1] - stepwidth, Tse0[2, 2], length.out=steps))
        
        GridB <- expand.grid(seq(Tse0[1, 2], Tse0[1, 2] + stepwidth, length.out=steps),
                             seq(Tse0[2, 1], Tse0[2, 2] + stepwidth, length.out=steps))
        
        GridC <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 2], length.out=steps),
                             seq(Tse0[2, 1] - stepwidth, Tse0[2, 1], length.out=steps))
        
        GridD <- expand.grid(seq(Tse0[1, 1], Tse0[1, 2] + stepwidth, length.out=steps),
                             seq(Tse0[2, 2], Tse0[2, 2] + stepwidth, length.out=steps))
        
        Grid <- rbind(GridA, GridB, GridC, GridD)
        
        FindcrTse <- apply(Grid, 1, function(x){
          theta <- x
          part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
        })
        
        CrTse <- cbind(Grid, FindcrTse)[FindcrTse==1, ]
        
        Tse <- t(apply(CrTse[, -(p + 1)], 2, range))
        
      }else{
        
        Tse <- Tse0
        
      }
      
    }
    
    TseOut <- cbind(est, Tse)
    rownames(TseOut) <- colnames(dat)
    colnames(TseOut) <- c("estimate", "lower", "upper")
    
    return(TseOut)
    
  }
  
  if(method=="tseng.brown"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:p){
      togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrTse <- apply(grid, 1, function(x){
      theta <- x
      (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
        qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
    })
    
    crTse <- cbind(grid, findcrTse)[findcrTse==1, ]
    
    if(nrow(crTse)==0){
      Tse0 <- matrix(rep(Inf, 2 * p), 2)
    }else{
      Tse0 <- t(apply(crTse[, -(p + 1)], 2, range))
    }
    
    if(Tse0[1, 1]==Inf | Tse0[1, 1]==-Inf){
      
      Tse <- matrix(rep(0, 2 * p), 2)
      
    }else{
      
      if(min(abs(Tse0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Tse0[, 2] - est - 2 * poolvar)) < 0.001){
        
        togrid2 <- list()
        
        for(i in 1:p){
          togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
        }
        
        grid2 <- expand.grid(togrid2)
        
        findcrTse2 <- apply(grid2, 1, function(x){
          theta <- x
          (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
            qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
        })
        
        crTse2 <- cbind(grid2, findcrTse2)[findcrTse2==1, ]
        
        if(min(crTse2[, 1])==min(grid2[, 1]) | max(crTse2[, 1])==max(grid2[, 1]) |
             min(crTse2[, 2])==min(grid2[, 2]) | max(crTse2[, 2])==max(grid2[, 2])){
          
          togrid3 <- list()
          
          for(i in 1:p){
            togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
          }
          
          grid3 <- expand.grid(togrid3)
          
          findcrTse3 <- apply(grid3, 1, function(x){
            theta <- x
            (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
              qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
          })
          
          crTse3 <- cbind(grid3, findcrTse3)[findcrTse3==1, ]
          
          crTse2 <- crTse3
          
        }
        
        Tse0 <- t(apply(crTse2[, -(p + 1)], 2, range))
        
      }
      
      stepwidth <- 4 * poolvar / steps
      
      if(p==2){
        
        GridA <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 1], length.out=steps),
                             seq(Tse0[2, 1] - stepwidth, Tse0[2, 2], length.out=steps))
        
        GridB <- expand.grid(seq(Tse0[1, 2], Tse0[1, 2] + stepwidth, length.out=steps),
                             seq(Tse0[2, 1], Tse0[2, 2] + stepwidth, length.out=steps))
        
        GridC <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 2], length.out=steps),
                             seq(Tse0[2, 1] - stepwidth, Tse0[2, 1], length.out=steps))
        
        GridD <- expand.grid(seq(Tse0[1, 1], Tse0[1, 2] + stepwidth, length.out=steps),
                             seq(Tse0[2, 2], Tse0[2, 2] + stepwidth, length.out=steps))
        
        Grid <- rbind(GridA, GridB, GridC, GridD)
        
        FindcrTse <- apply(Grid, 1, function(x){
          theta <- x
          (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
            qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
        })
        
        CrTse <- cbind(Grid, FindcrTse)[FindcrTse==1, ]
        
        Tse <- t(apply(CrTse[, -(p + 1)], 2, range))
        
      }else{
        
        Tse <- Tse0
        
      }
      
    }
    
    TseOut <- cbind(est, Tse)
    rownames(TseOut) <- colnames(dat)
    colnames(TseOut) <- c("estimate", "lower", "upper")
    
    return(TseOut)
    
  }
  
}
