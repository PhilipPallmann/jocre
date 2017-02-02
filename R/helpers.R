print.JOC <- summary.JOC <- function(x, digits=max(3, getOption("digits") - 4), ...){
  
  cat(paste("Parameter estimates and projected boundaries of the ", x$p, "-dimensional ", 100 * (1 - x$alpha),
            "% simultaneous confidence region\n", sep=""))
  
  res <- cbind(round(x$est, digits), round(x$ci, digits))
  rownames(res) <- colnames(x$dat)
  colnames(res) <- c("Estimate", "Lower", "Upper")
  
  print(res)
  
}

plot.JOC <- function(x, equi=log(c(0.8, 1.25)), axnames=NULL, main=NULL, xlim=log(c(0.77, 1.3)),
                     ylim=log(c(0.77, 1.3)), col="black", ...){
  
  if(nrow(x$ci)!=2){
    stop("Plotting only allowed for regions or intervals in 2 dimensions.")
  }
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(x$dat)
  }else{
    axisnames <- axnames
  }
  
  par(mar=c(5, 5, 4, 2))
  plot(0, xlim=xlim, ylim=ylim, las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main, ...)
  if(is.null(equi)==FALSE){
    if(length(equi)!=2){
      stop("Length of equi must be 2.")
    }
    rect(equi[1], equi[1], equi[2], equi[2], col="gray95", border=NA)
  }
  if(x$method %in% c("limacon.asy", "limacon.fin", "tseng", "tseng.brown")){
    points(x$cr, pch=20, col=col)
  }
  if(x$method %in% c("emp.bayes", "hotelling", "standard.cor", "standard.ind")){
    polygon(x$cr[chull(x$cr), -3], col=col, border=col)
  }
  if(x$method %in% c("expanded", "tost")){
    segments(x0=ci[1], x1=ci[3], y0=est[2], y1=est[2], lwd=2, col=col)
    segments(y0=ci[2], y1=ci[4], x0=est[1], x1=est[1], lwd=2, col=col)
  }
  points(est[1], est[2], pch=19, col="white")
  points(0, 0, pch="+", col="white", cex=2)
  par(mar=c(5, 4, 4, 2))
  
}

print.JOCMV <- summary.JOCMV <- function(x, digits=max(3, getOption("digits") - 4), ...){
  
  cat(paste("Parameter estimate and projected boundaries of the 2-dimensional ", 100 * (1 - x$alpha),
            "% simultaneous confidence region\n", sep=""))
  
  res <- cbind(round(x$est, digits), round(x$ci, digits))
  colnames(res) <- c("Estimate", "Lower", "Upper")
  
  print(res)
  
}

plot.JOCMV <- function(x, axnames=NULL, main=NULL, xlim=NULL, ylim=NULL, col="black", ...){
  
  if(is.null(xlim)==FALSE){
    xlims <- xlim
  }else{
    xlims <- range(x$cr[, 1])
  }
  
  if(is.null(ylim)==FALSE){
    ylims <- ylim
  }else{
    ylims <- range(x$cr[, 2])
  }
  
  if(is.null(axnames)==TRUE){
    if(x$scale=="var"){
      axisnames <- c("Mean", "Variance")
    }else{
      axisnames <- c("Mean", "SD")
    }
  }else{
    axisnames <- axnames
  }
  
  par(mar=c(5, 5, 4, 2))
  plot(0, xlim=xlims, ylim=ylims, las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
  polygon(x$cr[chull(x$cr[, ]), ], col=col, border=col)
  points(x$est, x$s^2, pch=19, col="white")
  par(mar=c(5, 4, 4, 2))
  
}