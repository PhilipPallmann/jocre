\name{jocre-package}

\alias{jocre-package}
\alias{jocre}

\docType{package}

\title{
Joint confidence regions
}

\description{
This package provides functions for computing and plotting joint confidence regions as well as (simultaneous) confidence intervals, with a focus on multivariate normal parameter vectors and parameters of the normal distribution i.e., mean and variance.}

\details{
\tabular{ll}{
Package: \tab jocre\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2016-09-27\cr
License: \tab GPL-2\cr
}
\code{plot2D} generates plots of joint confidence regions or (simultaneous) confidence intervals in two dimensions and allows for fine tuning e.g., shading of an equivalence region. The confidence region's boundaries can be found with \code{confset}, whereas \code{confints} computes (simultaneous) confidence intervals.

\code{plotMV2D} produces graphics of joint confidence regions for the mean and variance (or standard deviation) of a normal distribution.

Insights into the conservativeness of the two one-sided test procedure can be obtained with \code{iutsize}.

CIGE \dots
}
\author{
Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})
}
\references{
Philip Pallmann (2016) Joint confidence regions with the R package jocre. In preparation.
}

\keyword{package}

\examples{
data("marzo")
#plot2D()

data("wires")
#plotMV2D()
}
