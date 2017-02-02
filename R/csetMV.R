csetMV <- function(dat, n, method, alpha=0.1, scale="var", steps=400){
  
  method <- match.arg(method, choices=c("mood", "large", "plugin", "pluginF", "lrt", "cheng.iles", "min.area"))
  scale <- match.arg(scale, choices=c("var", "sd"))
  
  if((method %in% c("cheng.iles", "min.area")) & scale=="sd"){
    stop("Please choose scale='var' for this method.")
  }
  
  if(method %in% c("mood", "large", "plugin", "pluginF", "lrt")){
    
    if(is.vector(dat)!=TRUE){
      stop("dat must be a vector of numeric values.")
    }
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    searchwidth <- 2
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * n / qchisq(df=df, 1 - alpha/16),
                       s^2 * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    if(method=="mood"){
      
      grid[, 3] <- mea - qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) < grid[, 1] &
        grid[, 1] < mea + qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) &
        s^2 * n / qchisq(df=df, 1 - alpha/2) < grid[, 2] &
        grid[, 2] < s^2 * n / qchisq(df=df, alpha/2)
      
      crFinal <- grid[grid[, 3]==TRUE, ]
      
      while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
              min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
        
        togrid <- list()
        togrid[[1]] <- seq(mea - searchwidth * qnorm(1 - alpha/16) * s / sqrt(n),
                           mea + searchwidth * qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
        togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                           s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
        
        grid <- expand.grid(togrid)
        
        grid[, 3] <- mea - qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) < grid[, 1] &
          grid[, 1] < mea + qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) &
          s^2 * n / qchisq(df=df, 1 - alpha/2) < grid[, 2] &
          grid[, 2] < s^2 * n / qchisq(df=df, alpha/2)
        
        crFinal <- grid[grid[, 3]==TRUE, ]
        ciFinal <- t(apply(crFinal[, -3], 2, range, na.rm=TRUE))
        
        searchwidth <- searchwidth + 1
        
      }
      
    }
    
    if(method=="large"){
      
      grid[, 3] <- n/grid[, 2] * (mea - grid[, 1])^2 + n/(2 * grid[, 2]^2) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
      
      crFinal <- grid[grid[, 3]==TRUE, ]
      
      while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
              min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
        
        togrid <- list()
        togrid[[1]] <- seq(mea - searchwidth * qnorm(1 - alpha/16) * s / sqrt(n),
                           mea + searchwidth * qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
        togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                           s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
        
        grid <- expand.grid(togrid)
        
        grid[, 3] <- n/grid[, 2] * (mea - grid[, 1])^2 + n/(2 * grid[, 2]^2) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
        
        crFinal <- grid[grid[, 3]==TRUE, ]
        ciFinal <- t(apply(crFinal[, -3], 2, range, na.rm=TRUE))
        
        searchwidth <- searchwidth + 1
        
      }
      
    }
    
    if(method=="plugin"){
      
      grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
      
      crFinal <- grid[grid[, 3]==TRUE, ]
      
      while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
              min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
        
        togrid <- list()
        togrid[[1]] <- seq(mea - searchwidth * qnorm(1 - alpha/16) * s / sqrt(n),
                           mea + searchwidth * qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
        togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                           s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
        
        grid <- expand.grid(togrid)
        
        grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
        
        crFinal <- grid[grid[, 3]==TRUE, ]
        ciFinal <- t(apply(crFinal[, -3], 2, range, na.rm=TRUE))
        
        searchwidth <- searchwidth + 1
        
      }
      
    }
    
    if(method=="pluginF"){
      
      grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qf(1 - alpha, df1=2, df2=n - 2)
      
      crFinal <- grid[grid[, 3]==TRUE, ]
      
      while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
              min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
        
        togrid <- list()
        togrid[[1]] <- seq(mea - searchwidth * qnorm(1 - alpha/16) * s / sqrt(n),
                           mea + searchwidth * qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
        togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                           s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
        
        grid <- expand.grid(togrid)
        
        grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qf(1 - alpha, df1=2, df2=n - 2)
        
        crFinal <- grid[grid[, 3]==TRUE, ]
        ciFinal <- t(apply(crFinal[, -3], 2, range, na.rm=TRUE))
        
        searchwidth <- searchwidth + 1
        
      }
      
    }
    
    if(method=="lrt"){
      
      grid[, 3] <- n * log(grid[, 2] / s^2) + n * s^2 / grid[, 2] + n * (mea - grid[, 1])^2 / grid[, 2] - n < qchisq(1 - alpha, df=2)
      
      crFinal <- grid[grid[, 3]==TRUE, ]
      
      while(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
              min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
        
        togrid <- list()
        togrid[[1]] <- seq(mea - searchwidth * qnorm(1 - alpha/16) * s / sqrt(n),
                           mea + searchwidth * qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
        togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                           s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
        
        grid <- expand.grid(togrid)
        
        grid[, 3] <- n * log(grid[, 2] / s^2) + n * s^2 / grid[, 2] + n * (mea - grid[, 1])^2 / grid[, 2] - n < qchisq(1 - alpha, df=2)
        
        crFinal <- grid[grid[, 3]==TRUE, ]
        ciFinal <- t(apply(crFinal[, -3], 2, range, na.rm=TRUE))
        
        searchwidth <- searchwidth + 1
        
      }
      
    }
    
  }
  
  if(method=="cheng.iles"){
    
    int.cheng<-function(lob,hib,q,n,nint){
      
      bvec<-(1:nint-0.5)*(hib-lob)/nint+lob
      bdens<-double(nint)
      int<-double(nint)
      for (i in 1:nint){
        b<-bvec[i]
        a2<-sqrt((q*b*b-2*n*(1-b)*(1-b))/n)
        logbdens<-((n-1)/2)*log(n-1)-n*log(b)-(n-1)/(2*b*b)-lgamma((n-1)/2)-((n-3)/2)*log(2)
        bdens[i]<-exp(logbdens)
        aprob<-2*pnorm(a2,0,sd=b/sqrt(n))-1
        int[i]<-aprob*bdens[i]
      }
      cov.prob<-(sum(int))*(hib-lob)/nint
      
      cov.prob
      
    }
    
    int.chengr<-function(lob,hib,q,n,tol){
      
      oldbest<-0
      expo<-1
      nint<-5*(2**expo)
      best<-int.cheng(lob,hib,q,n,nint)
      err<-abs(best-oldbest)/(best+oldbest)
      while(err>tol){
        expo<-expo+1
        nint<-5*(2**expo)
        oldbest<-best
        best<-int.cheng(lob,hib,q,n,nint)
        err<-abs(best-oldbest)/(best+oldbest)
      }
      
      best
      
    }
    
    covprob.cheng<-function(n,q){
      
      tol<-0.00001
      
      a1<-2*n-q
      b1<--4*n
      c1<-2*n
      lob<-(-b1-sqrt(b1*b1-4*a1*c1))/(2*a1)
      hib<-(-b1+sqrt(b1*b1-4*a1*c1))/(2*a1)
      
      cp<-int.chengr(lob,hib,q,n,tol)
      
      list(cp=cp,lob=lob,hib=hib)
      
    }
    
    chengiles<-function(n,dcl){
      
      #Here we note a value that is too small.
      
      lo<-0
      loval<-0
      
      #Here we find a value that gives a coverage probability
      # higher than desired.
      
      hi<-1
      #hival<-covprob.cheng(n,hi)$cp
      
      #while (hival<dcl){
      #  hi<-hi*2
      #  hival<-covprob.cheng(n,hi)$cp}
      
      hi<-2*n
      hival<-1
      
      #Having found two bounding values, we now do bisection to
      # find the appropriate cutoff value.
      
      mid<-(lo+hi)/2
      midval<-covprob.cheng(n,mid)
      err<-abs(midval$cp-dcl)
      while (err>0.00001){
        if (midval$cp>dcl){
          hi<-mid
          hival<-midval$cp}
        if (midval$cp<dcl){
          lo<-mid
          loval<-midval$cp}
        mid<-(lo+hi)/2
        midval<-covprob.cheng(n,mid)
        err<-abs(midval$cp-dcl)
      }
      
      #The value we report is the cutoff.
      
      list(cutoff=mid,lob=midval$lob,hib=midval$hib)
      
    }
    
    plotcheng <-function(n,dcl,int){
      
      out<-chengiles(n,dcl)
      
      nint<-int
      pi<-4*atan(1)
      vec<-(sin(seq(-pi/2,pi/2,length=nint))+1)/2
      bvec<-out$lob+vec*(out$hib-out$lob)
      q<-out$cutoff
      a1<-double(nint)
      for (i in 2:(nint-1)){
        b<-bvec[i]
        a1[i]<--sqrt((q*b*b-2*n*(1-b)*(1-b))/n)
      }
      
      list(avec=a1,bvec=bvec)
      
    }
    
    forgrid <- plotcheng(n=n, dcl=1 - alpha, int=1000)
    crFinal <- cbind(c(forgrid$avec, -forgrid$avec), c(forgrid$bvec, forgrid$bvec))
    
  }
  
  if(method=="min.area"){
    
    covprob <- function(n, tar){
      
      tol<-0.00001
      bmax<-sqrt((n-1)/(n+1))
      
      #Here we find the limits in sigma.
      #First we look for the lower limit.
      
      lo<-0
      loval<--1
      hi<-bmax
      hival<-1
      mid<-(lo+hi)/2
      midval<--(n+1)*log(mid) - (n-1)/(2*mid*mid) - tar
      while (abs(hi-lo)>0.000001){
        if (midval>=0){
          hi<-mid
          hival<-midval}
        if (midval<0){
          lo<-mid
          loval<-midval}
        mid<-(lo+hi)/2
        midval<--(n+1)*log(mid) - (n-1)/(2*mid*mid) - tar
      }
      lob<-mid
      
      #Now we look for the upper limit.
      #First we find a value that exceeds the upper limit.
      
      lo<-bmax
      loval<-1
      hi<-bmax+1
      hival<--(n+1)*log(hi) - (n-1)/(2*hi*hi)-tar
      while (hival>0){
        hi<-hi*2
        hival<--(n+1)*log(hi) - (n-1)/(2*hi*hi)-tar}
      lo<-max(bmax,hi/2)
      mid<-(lo+hi)/2
      midval<--(n+1)*log(mid) - (n-1)/(2*mid*mid) - tar
      while (abs(hi-lo)>0.000001){
        if (midval>=0){
          lo<-mid
          loval<-midval}
        if (midval<0){
          hi<-mid
          hival<-midval}
        mid<-(lo+hi)/2
        midval<--(n+1)*log(mid) - (n-1)/(2*mid*mid) - tar}
      hib<-mid
      
      cp<-int.functr(lob,hib,tar,n,tol)
      list(cp=cp,lob=lob,hib=hib)
      
    }
    
    int.funct <- function(lob, hib, tar, n, nint){
      
      bvec<-(1:nint-0.5)*(hib-lob)/nint+lob
      a2<-double(nint)
      bdens<-double(nint)
      int<-double(nint)
      for (i in 1:nint){
        b<-bvec[i]
        a2[i]<-sqrt(2*b*b*(1/n)*(-tar-(n+1)*log(b)-(n-1)/(2*b*b)))
        logbdens<-((n-1)/2)*log(n-1)-n*log(b)-(n-1)/(2*b*b)-lgamma((n-1)/2)-((n-3)/2)*log(2)
        bdens[i]<-exp(logbdens)
        aprob<-2*pnorm(a2[i],0,sd=b/sqrt(n))-1
        int[i]<-aprob*bdens[i]
      }
      cov.prob<-(sum(int))*(hib-lob)/nint
      
      cov.prob
      
    }
    
    int.functr <- function(lob, hib, tar, n, tol){
      
      oldbest<-0
      expo<-1
      nint<-5*(2**expo)
      best<-int.funct(lob,hib,tar,n,nint)
      err<-abs(best-oldbest)/(best+oldbest)
      while(err>tol){
        expo<-expo+1
        nint<-5*(2**expo)
        oldbest<-best
        best<-int.funct(lob,hib,tar,n,nint)
        err<-abs(best-oldbest)/(best+oldbest)
      }
      best
      
    }
    
    minarea <- function(n, dcl){
      
      #Here we find the maximal value.
      
      bmax<-sqrt((n-1)/(n+1))
      maxval<--(n+1)*log(bmax) - (n-1)/(2*bmax*bmax)
      
      #Here we find a value that gives a coverage probability
      # smaller than desired.
      
      diff<-0
      cp<-0
      tol<-0.00001
      
      while (cp<dcl){
        
        diff<-diff+1
        tar<-maxval-diff
        
        cp<-covprob(n,tar)$cp
      }
      
      #Having found a large enough "diff" value, we now start
      # doing bisection to achieve the desired coverage level.
      
      hi<-maxval-diff+1
      hival<-0
      lo<-maxval-diff
      loval<-cp
      mid<-(lo+hi)/2
      midval<-covprob(n,mid)
      err<-abs(midval$cp-dcl)
      while (err>0.00001){
        if (midval$cp>dcl){
          lo<-mid
          loval<-midval$cp}
        if (midval$cp<dcl){
          hi<-mid
          hival<-midval$cp}
        mid<-(lo+hi)/2
        midval<-covprob(n,mid)
        err<-abs(midval$cp-dcl)
      }
      
      #The value we report is the cutoff.
      
      list(cutoff=mid,lob=midval$lob,hib=midval$hib)
      
    }
    
    plotci <- function(n, dcl, int){
      
      out<-minarea(n,dcl)
      
      nint<-int
      pi<-4*atan(1)
      vec<-(sin(seq(-pi/2,pi/2,length=nint))+1)/2
      bvec<-out$lob+vec*(out$hib-out$lob)
      a1<-double(nint)
      a2<-double(nint)
      for (i in 2:(nint-1)){
        b<-bvec[i]
        a1[i]<--sqrt(2*b*b*(1/n)*(-out$cutoff-(n+1)*log(b)-(n-1)/(2*b*b)))
        a2[i]<--a1[i]
      }
      
      list(avec=a1,bvec=bvec)
      
    }
    
    forgrid <- plotci(n=n, dcl=1 - alpha, int=1000)
    crFinal <- cbind(c(forgrid$avec, -forgrid$avec), c(forgrid$bvec, forgrid$bvec))
    
  }
  
  if(scale=="sd"){
    crFinal[, 2] <- sqrt(crFinal[, 2])
    ciFinal[2, ] <- sqrt(ciFinal[2, ])
  }
  
  Out <- list()
  
  Out$cr <- rbind(ddply(crFinal, .(Var1), summarise, Var2=range(Var2)),
                  ddply(crFinal, .(Var2), summarise, Var1=range(Var1)))
  Out$ci <- ciFinal
  Out$n <- n
  if(method %in% c("cheng.iles", "min.area")){
    Out$est <- 0
    Out$s <- 1
  }else{
    Out$est <- mea
    Out$s <- s
  }
  Out$df <- n - 1
  Out$scale <- scale
  Out$dat <- dat
  Out$method <- method
  Out$alpha <- alpha
  Out$steps <- steps
  
  class(Out) <- "JOCMV"
  
  return(Out)
  
}