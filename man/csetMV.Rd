\name{csetMV}

\alias{csetMV}

\title{Confidence regions and intervals for a normal mean and variance}

\description{Computes boundaries of (simultaneous) confidence regions for the mean and variance of a normal distribution using different methods.}

\usage{cset(dat, n, method, alpha=0.1, scale="var", steps=400)}

\arguments{
  \item{dat}{A vector of numeric values assumed to follow a normal distribution. Not required for \code{method="cheng.iles"} or \code{"min.area"}.}
  \item{n}{A numeric value giving the sample size. Only required for \code{method="cheng.iles"} or \code{"min.area"}.}
  \item{method}{A character string specifying the method to be used. See details for available \code{method}s.}
  \item{alpha}{A numeric value giving the type I error level to be controlled. Default is 0.1.}
  \item{scale}{A character string specifying whether the variance (\code{var}) or standard deviation (\code{sd}) is to be plotted on the y axis. Not required for \code{method="cheng.iles"} or \code{"min.area"}.}
  \item{steps}{An integer setting the initial number of steps for the search algorithm.}
}

\details{Available \code{method}s are: \code{mood} for the classical region described in Mood (1950); \code{large} for the large-sample approximation region described in section 4.1 of Arnold & Shavelle (1998); \code{plugin} for a plug-in variant of the large-sample approximation region described in section 4.2 of Arnold & Shavelle (1998); \code{pluginF} for the plug-in variant of the large-sample approximation region described in section 4.3 of Arnold & Shavelle (1998) using an asymptotic F distribution as in Douglas (1993); \code{lrt} for the likelihood ratio test region described in section 4.4 of Arnold & Shavelle (1998); \code{cheng.iles} for the region described in Cheng & Iles (1983); \code{min.area} for the minimum-area region described in Frey et al. (2009).}

\value{An object of class \code{JOCMV}.}

\references{
Barry C. Arnold & Robert M. Shaville (1998) Joint confidence sets for the mean and variance of a normal distribution. The American Statistician, 52(2), 133--140.

R. C. H. Cheng & T. C. Iles (1983) Confidence bands for cumulative distribution functions of continuous random variables. Technometrics, 25(1), 77--86.

J. B. Douglas (1993) Confidence regions for parameter pairs. The American Statistician, 47(1), 43--45.

Jesse Frey, Osvaldo Marrero, Douglas Norton (2009) Minimum-area confidence sets for a normal distribution. Journal of Statistical Planning and Inference, 139(3), 1023--1032.

Alexander M. Mood (1950) Introduction to the Theory of Statistics. McGraw-Hill, New York, NY.
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

\note{Warning: please use with care! Some of the functionality has not yet been thoroughly tested.}

\seealso{\code{\link{cset}} for (simultaneous) confidence regions and intervals around multivariate normal means.}

\examples{
\dontrun{
# Simultaneous 90% confidence regions for univariate normal data

univar <- rnorm(n=50)

moodvar <- csetMV(dat=univar, method="mood", alpha=0.1, scale="var")
moodsd <- csetMV(dat=univar, method="mood", alpha=0.1, scale="sd")

plot(moodvar)
plot(moodsd)
}
}

\keyword{multivariate}