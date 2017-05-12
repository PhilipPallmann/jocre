cset <- function(dat, method, alpha=0.1, steps=NULL, nboot=1e4, TsengBrownA=1, TsengBrownB=1){
  
  Var1 <- Var2 <- NULL # just to appease RCMD check
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  est <- matrix(colMeans(dat), p)
  sd <- matrix(apply(dat, 2, sd), p)
  poolvar <- var(as.vector(as.matrix(dat)))
  cov <- cov(dat)
  solved <- solve(cov)
  #s2 <- poolvar/n
  #s <- sqrt(poolvar/n)
  
  method <- match.arg(method, choices=c("boot.kern", "emp.bayes", "expanded", "fix.seq", "hotelling",
                                        "limacon.asy", "limacon.fin", "standard.cor", "standard.ind",
                                        "tost", "tseng", "tseng.brown"))
  
  if(is.null(steps)==TRUE){
    if(p==2){
      steps <- 300
    }
    if(p > 2){
      steps <- 50
    }
  }
  
  if(method=="boot.kern"){
    
    if(p > 2){stop("The fixed sequence procedure is currently only implemented for 2-dimensional data.")}
    
    bivarmean <- function(x, d) {
      e <- x[d, ]
      return(c(mean(e[, 1]), mean(e[, 2])))
    } # werden 1. und 2. Spalte immer zusammen gesampelt?!
    
    b <- boot(dat, bivarmean, R=nboot)
    bdat <- as.data.frame(b$t)
    
    kern <- bkde2D(bdat, bandwidth=sapply(bdat, dpik))
    
    alphasum <- sum(kern$fhat) * alpha
    
    while(sum(kern$fhat, na.rm=TRUE) > alphasum){
      kern$fhat[which.max(kern$fhat)] <- NA
    }
    
    KERN <- is.na(kern$fhat)
    K <- expand.grid(kern$x1, kern$x2)
    K$truefalse <- as.vector(KERN)
    K <- K[K$truefalse==1, ]
    
    ciFinal <- rbind(range(K$Var1), range(K$Var2))
    
    hu <- chull(K[, -3])
    hull <- c(hu, hu[1])
    
    crFinal <- K[hull, -3]
    
    #par(mar=c(5, 5, 4, 2))
    #plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
    #     cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    #if(is.null(equi)==FALSE){
    #  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    #}
    #polygon(K[hull, ], col=col)
    #points(est[1], est[2], pch=19, col="white")
    #par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method=="emp.bayes"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
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
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        sqrt(sum((theta - JSplus)^2)) < (sqrt(poolvar/n) * sqrt(vE2))
      })
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- matrix(x, p)
        sqrt(sum((theta - JSplus)^2)) < (sqrt(poolvar/n) * sqrt(vE2))
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- matrix(x, p)
        sqrt(sum((theta - JSplus)^2)) < (sqrt(poolvar/n) * sqrt(vE2))
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="expanded"){
    
    ciFinal <- matrix(c(est - sd * qt(1 - alpha, df) / sqrt(n),
                        est + sd * qt(1 - alpha, df) / sqrt(n)), ncol(dat))
    
    ciFinal[, 1] <- ifelse(ciFinal[, 1] > 0, 0, ciFinal[, 1])
    ciFinal[, 2] <- ifelse(ciFinal[, 2] < 0, 0, ciFinal[, 2])
    
    crFinal <- NULL
    
  }
  
  if(method=="fix.seq"){
    
    if(p > 2){
      stop("The fixed sequence procedure is currently only implemented for 2-dimensional data.")
    }
    
    t_1st <- c(min(0, est[1] - sqrt(cov[1, 1]) * qt(1 - alpha, df) / sqrt(n)),
               max(0, est[1] + sqrt(cov[1, 1]) * qt(1 - alpha, df) / sqrt(n)))
    
    t_2nd <- c(min(0, est[2] - sqrt(cov[2, 2]) * qt(1 - alpha, df) / sqrt(n)),
               max(0, est[2] + sqrt(cov[2, 2]) * qt(1 - alpha, df) / sqrt(n)))
    
    if(t_1st[1] > log(0.8) & t_1st[2] < log(1.25)){
      
      if(t_2nd[1] > log(0.8) & t_2nd[2] < log(1.25)){
        
        # both bioequivalent
        
        ma <- max(c(max(abs(t_1st)), max(abs(t_2nd))))
        T_1st <- T_2nd <- c(-ma, ma)
        
      }else{
        
        # only 1st bioequivalent
        
        T_1st <- c(log(0.8), log(1.25))
        T_2nd <- c(min(log(0.8), t_2nd[1]), max(log(1.25), t_2nd[2]))
        
      }
      
    }else{
      
      # 1st not bioequivalent (thus 2nd not tested)
      
      T_1st <- c(min(log(0.8), t_1st[1]), max(log(1.25), t_1st[2]))
      T_2nd <- c(NA, NA)
      
    }
    
    ciFinal <- rbind(T_1st, T_2nd)
    
    crFinal <- NULL
    
  }
  
  if(method=="hotelling"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    rhs <- qf(p=1 - alpha, df1=p, df2=df - p + 1) * p * df / (df - p + 1)
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- matrix(x, p)
        (n * t(est - theta) %*% solved %*% (est - theta)) < rhs 
      })
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- matrix(x, p)
        (n * t(est - theta) %*% solved %*% (est - theta)) < rhs 
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- matrix(x, p)
        (n * t(est - theta) %*% solved %*% (est - theta)) < rhs 
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="limacon.asy"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
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
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% solved %*% est) / sqrt(t(theta) %*% solved %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solved %*% theta) * sqrt(n)) 
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% solved %*% est) / sqrt(t(theta) %*% solved %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solved %*% theta) * sqrt(n)) 
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="limacon.fin"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
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
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="standard.cor"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    rhs <- (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=df))
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        t(est - theta) %*% solved %*% (est - theta) < rhs
      })
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- x
        t(est - theta) %*% solved %*% (est - theta) < rhs
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- x
        t(est - theta) %*% solved %*% (est - theta) < rhs
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="standard.ind"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    rhs <- (poolvar/n * p * qf(p=1 - alpha, df1=p, df2=df))
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
      togrid <- list()
      
      for(i in 1:p){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcr <- apply(grid, 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < rhs
      })
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
    }
    
    a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
    
    stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
    
    b <- a
    b[, 1] <- b[, 1] - stepwidth
    b[, 2] <- b[, 2] + stepwidth
    
    togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
    crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
    
    for(i in 1:p){
      
      togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
      togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
      
      if(p==2){
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i==2){
          
          biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
          
        }
        
      }else{
        
        biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
        biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
        
        if(i > 1 & i < p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
        }
        
        if(i==p){
          biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
          biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
        }
        
      }
      
      colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
      biglistA[[i]] <- as.matrix(biglistA[[i]])
      biglistB[[i]] <- as.matrix(biglistB[[i]])
      
      findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < rhs
      })
      
      findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < rhs
      })
      
      crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
      crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
      
      if(is.matrix(crFinalA[[i]])==FALSE){
        crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
      }
      
      if(is.matrix(crFinalB[[i]])==FALSE){
        crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
      }
      
      ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
      ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
      
    }
    
    ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
    
    FinalA <- ldply(crFinalA, data.frame)
    FinalB <- ldply(crFinalB, data.frame)
    
    colnames(crFinalX) <- colnames(FinalA)
    
    crFinal <- rbind(crFinalX, FinalA, FinalB)
    
  }
  
  if(method=="tost"){
    
    ciFinal <- matrix(c(est - sd * qt(1 - alpha, df) / sqrt(n),
                        est + sd * qt(1 - alpha, df) / sqrt(n)), ncol(dat))
    
    crFinal <- NULL
    
  }
  
  if(method=="tseng"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
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
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
      problem <- FALSE
      
      if(nrow(crFinalX)==0){
        ciFinal <- matrix(rep(Inf, 2 * p), 2)
        crFinal <- crFinalX
        problem <- TRUE
        break
      }
      
      if(ciFinalX[1, 1]==Inf | ciFinalX[1, 1]==-Inf){
        ciFinal <- matrix(rep(0, 2 * p), 2)
        crFinal <- crFinalX
        problem <- TRUE
        break
      }
      
    }
    
    if(problem!=TRUE){
      
      a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
      
      b <- a
      b[, 1] <- b[, 1] - stepwidth
      b[, 2] <- b[, 2] + stepwidth
      
      togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
      crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
      
      for(i in 1:p){
        
        togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
        togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
        
        if(p==2){
          
          biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
          biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
          
          if(i==2){
            
            biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
            
          }
          
        }else{
          
          biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
          biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
          
          if(i > 1 & i < p){
            biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
          }
          
          if(i==p){
            biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
          }
          
        }
        
        colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
        biglistA[[i]] <- as.matrix(biglistA[[i]])
        biglistB[[i]] <- as.matrix(biglistB[[i]])
        
        findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
          theta <- x
          lhs > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / (poolvar/n)))
        })
        
        findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
          theta <- x
          lhs > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / (poolvar/n)))
        })
        
        crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
        crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
        
        if(is.matrix(crFinalA[[i]])==FALSE){
          crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
        }
        
        if(is.matrix(crFinalB[[i]])==FALSE){
          crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
        }
        
        ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
        ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
        
      }
      
      ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
      
      FinalA <- ldply(crFinalA, data.frame)
      FinalB <- ldply(crFinalB, data.frame)
      
      colnames(crFinalX) <- colnames(FinalA)
      
      crFinal <- rbind(crFinalX, FinalA, FinalB)
      
    }
    
  }
  
  if(method=="tseng.brown"){
    
    searchwidth <- 1
    ciFinalX <- cbind(est - searchwidth/2 * poolvar, est + searchwidth/2 * poolvar)
    
    while(min(abs(ciFinalX[, 1] - (est - searchwidth/2 * poolvar))) < 0.001 | 
            min(abs(ciFinalX[, 2] - (est + searchwidth/2 * poolvar))) < 0.001){
      
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
      
      crFinalX <- cbind(grid, findcr)[findcr==1, ]
      ciFinalX <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      searchwidth <- 2 * searchwidth
      
      problem <- FALSE
      
      if(nrow(crFinalX)==0){
        ciFinal <- matrix(rep(Inf, 2 * p), 2)
        crFinal <- crFinalX
        problem <- TRUE
        break
      }
      
      if(ciFinalX[1, 1]==Inf | ciFinalX[1, 1]==-Inf){
        ciFinal <- matrix(rep(0, 2 * p), 2)
        crFinal <- crFinalX
        problem <- TRUE
        break
      }
      
    }
    
    if(problem!=TRUE){
      
      a <- t(apply(crFinalX[, -(p + 1)], 2, range, na.rm=TRUE))
      
      stepwidth <- diff(range(grid[, 1]))/steps # same for all variables
      
      b <- a
      b[, 1] <- b[, 1] - stepwidth
      b[, 2] <- b[, 2] + stepwidth
      
      togridA <- togridB <- biglistA <- biglistB <- findcrA <- findcrB <- list()
      crFinalA <- crFinalB <- ciFinalA <- ciFinalB <- list()
      
      for(i in 1:p){
        
        togridA[[i]] <- seq(a[i, 1], b[i, 1], length.out=steps)
        togridB[[i]] <- seq(a[i, 2], b[i, 2], length.out=steps)
        
        if(p==2){
          
          biglistA[[i]] <- expand.grid(c(list(togridA[[i]], crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
          biglistB[[i]] <- expand.grid(c(list(togridB[[i]], crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
          
          if(i==2){
            
            biglistA[[i]] <- cbind(biglistA[[i]][, 2], biglistA[[i]][, 1])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2], biglistB[[i]][, 1])
            
          }
          
        }else{
          
          biglistA[[i]] <- expand.grid(c(list(togridA[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 1], ][, -c(i, p + 1)])))
          biglistB[[i]] <- expand.grid(c(list(togridB[[i]]), as.list(crFinalX[crFinalX[, i]==a[i, 2], ][, -c(i, p + 1)])))
          
          if(i > 1 & i < p){
            biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1], biglistA[[i]][, (i + 1):p])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1], biglistB[[i]][, (i + 1):p])
          }
          
          if(i==p){
            biglistA[[i]] <- cbind(biglistA[[i]][, 2:i], biglistA[[i]][, 1])
            biglistB[[i]] <- cbind(biglistB[[i]][, 2:i], biglistB[[i]][, 1])
          }
          
        }
        
        colnames(biglistA[[i]]) <- colnames(biglistB[[i]]) <- 1:p
        biglistA[[i]] <- as.matrix(biglistA[[i]])
        biglistB[[i]] <- as.matrix(biglistB[[i]])
        
        findcrA[[i]] <- apply(biglistA[[i]], 1, function(x){
          theta <- x
          (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
            qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
        })
        
        findcrB[[i]] <- apply(biglistB[[i]], 1, function(x){
          theta <- x
          (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
            qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
        })
        
        crFinalA[[i]] <- cbind(biglistA[[i]], findcrA[[i]])[findcrA[[i]]==1, ]
        crFinalB[[i]] <- cbind(biglistB[[i]], findcrB[[i]])[findcrB[[i]]==1, ]
        
        if(is.matrix(crFinalA[[i]])==FALSE){
          crFinalA[[i]] <- rbind(crFinalA[[i]], crFinalA[[i]])
        }
        
        if(is.matrix(crFinalB[[i]])==FALSE){
          crFinalB[[i]] <- rbind(crFinalB[[i]], crFinalB[[i]])
        }
        
        ciFinalA[[i]] <- t(apply(crFinalA[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 1]
        ciFinalB[[i]] <- t(apply(crFinalB[[i]][, -(p + 1)], 2, range, na.rm=TRUE))[i, 2]
        
      }
      
      ciFinal <- matrix(unlist(c(ciFinalA, ciFinalB)), nrow=p)
      
      FinalA <- ldply(crFinalA, data.frame)
      FinalB <- ldply(crFinalB, data.frame)
      
      colnames(crFinalX) <- colnames(FinalA)
      
      crFinal <- rbind(crFinalX, FinalA, FinalB)
      
    }
    
  }
  
  Out <- list()
  
  if(p==2 & is.null(crFinal)==FALSE){
    colnames(crFinal) <- c("Var1", "Var2")
    Out$cr <- rbind(ddply(crFinal, .(Var1), summarise, Var2=range(Var2)),
                    ddply(crFinal, .(Var2), summarise, Var1=range(Var1)))
  }
  if(p > 2){
    Out$cr <- NULL
  }
  Out$ci <- ciFinal
  Out$n <- n
  Out$p <- p
  if(method=="emp.bayes"){
    Out$est <- JSplus
  }else{
    Out$est <- est
  }
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
