plot2D <- function(dat, method, alpha=0.1, equi=log(1.25), axnames=NULL, main="Title",
                   xlim=log(c(0.77, 1.3)), ylim=log(c(0.77, 1.3)), col="black", steps=400,
                   nboot=1e4, TsengBrownA=1, TsengBrownB=1){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  dat <- as.data.frame(dat)
  
  method <- match.arg(method, choices=c("bootkern", "emp.bayes", "expanded", "fixseq", "hotelling", "limacon.asy",
                                        "limacon.fin", "standard.cor", "standard.ind", "tost", "tseng", "tseng.brown"))
  
  searchwidth <- 2
  
  if(method=="bootkern"){
    
    est <- matrix(colMeans(dat), 2)
    
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
    
    cint <- rbind(range(K$Var1), range(K$Var2))
    
    hu <- chull(K[, -3])
    hull <- c(hu, hu[1])
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    polygon(K[hull, ], col=col)
    points(est[1], est[2], pch=19, col="white")
    par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method=="emp.bayes"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    s <- sqrt(poolvar / n)
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrCas <- apply(grid, 1, function(x){
      theta <- matrix(x, 2)
      sqrt(sum((theta - est)^2)) < (s * sqrt(2 * qf(p=1 - alpha, df1=2, df2=df)))
    })
    
    crFinal <- cbind(grid, findcrCas)[findcrCas==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrCas <- apply(grid, 1, function(x){
        theta <- matrix(x, 2)
        sqrt(sum((theta - est)^2)) < (s * sqrt(2 * qf(p=1 - alpha, df1=2, df2=df)))
      })
      
      crFinal <- cbind(grid, findcrCas)[findcrCas==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="expanded"){
    
    n <- nrow(dat)
    
    est <- matrix(colMeans(dat), 2)
    sd <- c(sd(dat[, 1]), sd(dat[, 2]))
    
    ci <- matrix(c(est - sd * qt(1 - alpha, n - 1) / sqrt(n),
                   est + sd * qt(1 - alpha, n - 1) / sqrt(n)), ncol(dat))
    ci[, 1] <- ifelse(ci[, 1] > 0, 0, ci[, 1])
    ci[, 2] <- ifelse(ci[, 2] < 0, 0, ci[, 2])
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    segments(x0=ci[1], x1=ci[3], y0=est[2], y1=est[2], lwd=2, col=col)
    segments(y0=ci[2], y1=ci[4], x0=est[1], x1=est[1], lwd=2, col=col)
    par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method=="fixseq"){
    
    n <- nrow(dat)
    p <- ncol(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), p)
    cov <- cov(dat)
    
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
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    segments(x0=T_1st[1], x1=T_1st[2], y0=est[2], y1=est[2], lwd=2, col=col)
    segments(y0=T_2nd[1], y1=T_2nd[2], x0=est[1], x1=est[1], lwd=2, col=col)
    par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method=="hotelling"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrHot <- apply(grid, 1, function(x){
      theta <- matrix(x, 2)
      (n * t(est - theta) %*% solve(cov) %*% (est - theta)) < qf(p=1 - alpha, df1=2, df2=df - 1) * 2 * df / (df - 1)
    })
    
    crFinal <- cbind(grid, findcrHot)[findcrHot==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrHot <- apply(grid, 1, function(x){
        theta <- matrix(x, 2)
        (n * t(est - theta) %*% solve(cov) %*% (est - theta)) < qf(p=1 - alpha, df1=2, df2=df - 1) * 2 * df / (df - 1)
      })
      
      crFinal <- cbind(grid, findcrHot)[findcrHot==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="limacon.asy"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrLim <- apply(grid, 1, function(x){
      theta <- matrix(x, 2)
      ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
        (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
    })
    
    crFinal <- cbind(grid, findcrLim)[findcrLim==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrLim <- apply(grid, 1, function(x){
        theta <- matrix(x, 2)
        ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
      })
      
      crFinal <- cbind(grid, findcrLim)[findcrLim==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }

  if(method=="limacon.fin"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    cov <- cov(dat)
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrLim <- apply(grid, 1, function(x){
      theta <- matrix(x, 2)
      ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
        ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
    })
    
    crFinal <- cbind(grid, findcrLim)[findcrLim==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrLim <- apply(grid, 1, function(x){
        theta <- matrix(x, 2)
        ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
          ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
      })
      
      crFinal <- cbind(grid, findcrLim)[findcrLim==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="standard.cor"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    LambdaInv <- solve(cov(dat))
    poolvar <- var(as.vector(as.matrix(dat)))
    #s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrUsu <- apply(grid, 1, function(x){
      theta <- x
      #sqrt(sum((est - theta)^2))^2 < (s2 * 2 * qf(p=1 - alpha, df1=2, df2=n - 1))
      t(est - theta) %*% LambdaInv %*% (est - theta) < (2 / (n - 1) * qf(p=1 - alpha, df1=2, df2=n - 1))
    })
    
    crFinal <- cbind(grid, findcrUsu)[findcrUsu==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrUsu <- apply(grid, 1, function(x){
        theta <- x
        #sqrt(sum((est - theta)^2))^2 < (s2 * 2 * qf(p=1 - alpha, df1=2, df2=n - 1))
        t(est - theta) %*% LambdaInv %*% (est - theta) < (2 / (n - 1) * qf(p=1 - alpha, df1=2, df2=n - 1))
      })
      
      crFinal <- cbind(grid, findcrUsu)[findcrUsu==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="standard.ind"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrUsu <- apply(grid, 1, function(x){
      theta <- x
      sqrt(sum((est - theta)^2))^2 < (s2 * 2 * qf(p=1 - alpha, df1=2, df2=n - 1))
    })
    
    crFinal <- cbind(grid, findcrUsu)[findcrUsu==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrUsu <- apply(grid, 1, function(x){
        theta <- x
        sqrt(sum((est - theta)^2))^2 < (s2 * 2 * qf(p=1 - alpha, df1=2, df2=n - 1))
      })
      
      crFinal <- cbind(grid, findcrUsu)[findcrUsu==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="tost"){
    
    n <- nrow(dat)
    
    est <- matrix(colMeans(dat), 2)
    sd <- c(sd(dat[, 1]), sd(dat[, 2]))
    
    ci <- c(est - sd * qt(1 - alpha, n - 1) / sqrt(n),
            est + sd * qt(1 - alpha, n - 1) / sqrt(n))
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    segments(x0=ci[1], x1=ci[3], y0=est[2], y1=est[2], lwd=2, col=col)
    segments(y0=ci[2], y1=ci[4], x0=est[1], x1=est[1], lwd=2, col=col)
    par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method=="tseng"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrTse <- apply(grid, 1, function(x){
      theta <- x
      (sqrt(sum(est^2))^2 / (2 * s2)) > qf(p=alpha, df1=2, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
    })
    
    crFinal <- cbind(grid, findcrTse)[findcrTse==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrTse <- apply(grid, 1, function(x){
        theta <- x
        (sqrt(sum(est^2))^2 / (2 * s2)) > qf(p=alpha, df1=2, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
      })
      
      crFinal <- cbind(grid, findcrTse)[findcrTse==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method=="tseng.brown"){
    
    n <- nrow(dat)
    df <- n - 1
    
    est <- matrix(colMeans(dat), 2)
    poolvar <- var(as.vector(as.matrix(dat)))
    s2 <- poolvar / n
    
    togrid <- list()
    
    for(i in 1:2){
      togrid[[i]] <- seq(est[i] - poolvar, est[i] + poolvar, length.out=steps)
    }
    
    grid <- expand.grid(togrid)
    
    findcrTse <- apply(grid, 1, function(x){
      theta <- x
      (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
        qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
    })
    
    crFinal <- cbind(grid, findcrTse)[findcrTse==1, ]
    
    while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
            min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      
      togrid <- list()
      
      for(i in 1:2){
        togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
      }
      
      grid <- expand.grid(togrid)
      
      findcrTse <- apply(grid, 1, function(x){
        theta <- x
        (sqrt(sum((est - theta * (1 + (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))))^2))^2) <
          qchisq(p=alpha, df=2, ncp=((sqrt(sum(theta^2))^2) * (1/(TsengBrownA + TsengBrownB * (sqrt(sum(theta^2))^2)))^2))
      })
      
      crFinal <- cbind(grid, findcrTse)[findcrTse==1, ]
      
      searchwidth <- 2 * searchwidth
      
    }
    
  }
  
  if(method %in% c("emp.bayes", "hotelling", "standard.cor", "standard.ind")){
    
    if(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
         min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      warning("The search grid is too narrow, please increase searchwidth.")
    }
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    polygon(crFinal[chull(crFinal[, -3]), -3], col=col, border=col)
    points(est[1], est[2], pch=19, col="white")
    par(mar=c(5, 4, 4, 2))
    
  }
  
  if(method %in% c("limacon.asy", "limacon.fin", "tseng", "tseng.brown")){
    
    if(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
         min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
      warning("The search grid is too narrow, please increase searchwidth.")
    }
    
    if(is.null(axnames)==TRUE){
      axisnames <- colnames(dat)
    }else{
      axisnames <- axnames
    }
    
    par(mar=c(5, 5, 4, 2))
    plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
         cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
    if(is.null(equi)==FALSE){
      rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
    }
    points(crFinal[, -3], pch=20, col=col)
    points(est[1], est[2], pch=19, col="white")
    par(mar=c(5, 4, 4, 2))
    
  }
  
}
