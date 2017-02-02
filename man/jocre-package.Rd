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
Version: \tab 0.3\cr
Date: \tab 2017-02-02\cr
License: \tab GPL-2\cr
}
\code{cset} computes joint confidence regions and (simultaneous) confidence intervals around multivariate normal means. The parameter estimates and interval bounds (in the case of confidence regions their boundaries are projected onto the axes) are displayed with \code{print} or \code{summary}, and \code{plot} generates a graph of the estimate and region or intervals and allows for some fine tuning e.g., shading of an equivalence region. 

\code{csetMV} computes joint confidence regions for the mean and variance (or standard deviation) of a normal distribution. \code{plot}, \code{print}, and \code{summary} produce the corresponding graphics and summaries.

Insights into the conservativeness of the two one-sided test procedure can be obtained with \code{iutsize}.

Confidence levels of a joint ellipse and marginal intervals can be "translated" back and forth with \code{translate}.
}
\author{
Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})
}
\references{
Philip Pallmann (2017) Joint confidence regions with the R package jocre. In preparation.
}

\keyword{package}

\examples{
data("marzo")
#plot(cset(dat=marzo, method="limacon"))

data("wires")
#plot(csetMV(dat=wires, method="mood"))
}
