\name{plotMV2D}

\alias{plotMV2D}

\title{Plotting 2D confidence sets for normal mean and variance}

\description{Creates graphs of (simultaneous) confidence regions for the mean and variance of a normal distribution using different methods.}

\usage{plotMV2D(dat, n, method, alpha=0.1, scale="var", axnames=NULL,
                main="Title", xlim=NULL, ylim=NULL, col="black", steps=400)}

\arguments{
  \item{dat}{A vector of numeric values assumed to follow a normal distribution. Not required for \code{method="cheng.iles"} or "min.area".}
  \item{n}{A numeric value giving the sample size. Only required for \code{method="cheng.iles"} or "min.area".}
  \item{method}{A character string specifying the method to be used. See details for available \code{method}s.}
  \item{alpha}{A numeric value giving the type I error level to be controlled. Default is \code{0.1}.}
  \item{scale}{A character string specifying whether the variance (\code{var}) or standard deviation (\code{sd}) is to be plotted on the y axis.}
  \item{axnames}{A vector of two character strings giving the x and y axis labels. Default is \code{NULL}.}
  \item{main}{A character string giving the plot title.}
  \item{xlim}{A numeric vector of length two specifying the plotting range on the x-axis. Default is \code{NULL}.}
  \item{ylim}{A numeric vector of length two specifying the plotting range on the y-axis. Default is \code{NULL}.}
  \item{col}{A character string specifying the colour of the plotted region or intervals.}
  \item{steps}{An integer setting the initial number of steps for the search algorithm.}
}

\details{xxx

Available \code{method}s are: \code{mood} for the classical region described in Mood (1950); \code{large} for the large-sample approximation region described in section 4.1 of Arnold & Shavelle (1998); \code{plugin} for a plug-in variant of the large-sample approximation region described in section 4.2 of Arnold & Shavelle (1998); \code{pluginF} for the plug-in variant of the large-sample approximation region described in section 4.3 of Arnold & Shavelle (1998) using an asymptotic F distribution as in Douglas (1993); \code{lrt} for the likelihood ratio test region described in section 4.4 of Arnold & Shavelle (1998); \code{cheng.iles} for the xxxxx described in Cheng & Iles (1983); \code{min.area} for the minimum-area region described in Frey et al. (2009).
}

\value{A graphical display of a (simultaneous) confidence set in two dimensions.}

\references{
Barry C. Arnold & Robert M. Shaville (1998) Joint confidence sets for the mean and variance of a normal distribution. The American Statistician, 52(2), 133--140.

R. C. H. Cheng & T. C. Iles (1983) Confidence bands for cumulative distribution functions of continuous random variables. Technometrics, 25(1), 77--86.

J. B. Douglas (1993) Confidence regions for parameter pairs. The American Statistician, 47(1), 43--45.

Jesse Frey, Osvaldo Marrero, Douglas Norton (2009) Minimum-area confidence sets for a normal distribution. Journal of Statistical Planning and Inference, 139(3), 1023--1032.

Alexander M. Mood (1950) Introduction to the Theory of Statistics. McGraw-Hill, New York, NY.
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk}) using code from xxxxx}

\note{x}

\seealso{x % objects to See Also as \code{\link{help}}
}

\examples{
# Generate normal data

mydata <- rnorm(n=50)

# Plot simultaneous 90% confidence regions

#plotMV2D(dat=mydata, method="mood", alpha=0.1, main="Mood")
#plotMV2D(dat=mydata, method="large", alpha=0.1, main="Large-sample")
#plotMV2D(dat=mydata, method="plugin", alpha=0.1, main="Plug-in")
#plotMV2D(dat=mydata, method="pluginF", alpha=0.1, main="Plug-in (F distribution)")
#plotMV2D(dat=mydata, method="lrt", alpha=0.1, main="Likelihood ratio test")

# more examples with cheng.iles and min.area
}

\keyword{multivariate}
