\name{helpers}

\alias{plot.JOC}
\alias{print.JOC}
\alias{summary.JOC}
\alias{plot.JOCMV}
\alias{print.JOCMV}
\alias{summary.JOCMV}

\title{Different generic functions for classes \code{JOC} and \code{JOCMV}}

\description{Generic functions for summarising and plotting objects of class \code{JOC} or \code{JOCMV}.}

\usage{
\method{plot}{JOC}(x, equi=log(c(0.8, 1.25)), axnames=NULL, main=NULL, xlim=log(c(0.77, 1.3)),
     ylim=log(c(0.77, 1.3)), col="black", \dots)

\method{print}{JOC}(x, digits=max(3, getOption("digits") - 4), \dots)

\method{summary}{JOC}(object, digits=max(3, getOption("digits") - 4), \dots)

\method{plot}{JOCMV}(x, axnames=NULL, main=NULL, xlim=NULL, ylim=NULL, col="black", \dots)

\method{print}{JOCMV}(x, digits=max(3, getOption("digits") - 4), \dots)

\method{summary}{JOCMV}(object, digits=max(3, getOption("digits") - 4), \dots)
}

\arguments{
  \item{x}{An output object of class \code{JOC} or \code{JOCMV}.}
  \item{object}{An output object of class \code{JOC} or \code{JOCMV}.}
  \item{digits}{A numeric value giving the number of significant digits to be printed.}
  \item{equi}{A numeric vector of length 2 specifying the equivalence region (lower and upper equivalence threshold) to be shaded in grey. When set to \code{NULL} no equivalence region is drawn. Default is \code{log(c(0.8, 1.25))}.}
  \item{axnames}{A vector of two character strings giving the x and y axis labels. For \code{plot.JOC}, when set to \code{NULL} the column names of \code{dat} are used as axis labels. Default is \code{NULL}.}
  \item{main}{A character string giving the plot title. Default is \code{NULL}.}
  \item{xlim}{A numeric vector of length two specifying the plotting range on the x-axis. Default is \code{log(c(0.77, 1.3))} for \code{plot.JOC} and \code{NULL} for \code{plot.JOCMV}.}
  \item{ylim}{A numeric vector of length two specifying the plotting range on the y-axis. Default is \code{log(c(0.77, 1.3))} for \code{plot.JOC} and \code{NULL} for \code{plot.JOCMV}.}
  \item{col}{A character string specifying the colour of the plotted region or intervals.}
  \item{\dots}{Further plotting arguments to be passed to methods. Type \code{?plot} for details.}
}

\details{\code{print} and \code{summary} summarise the estimates and confidence set boundaries of an object of class \code{JOC} or \code{JOCMV} that was created with \code{cset} or \code{csetMV}, respectively. \code{plot} displays a (simultaneous) confidence region or intervals when applied to an object of class \code{JOC} or \code{JOCMV}.}

\value{An on-screen summary or graphical display.}

\references{
Philip Pallmann & Thomas Jaki (2017) Simultaneous confidence regions and intervals for multivariate bioequivalence. Submitted to Statistics in Medicine.
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

\note{Warning: please use with care! Some of the functionality has not yet been thoroughly tested.}

\seealso{\code{\link{cset}} and \code{\link{csetMV}} for computing (simultaneous) confidence regions and intervals.}

\examples{
\dontrun{
# Example 1: simultaneous 90\% confidence region for bivariate data

bivar <- mvtnorm::rmvnorm(n=100, mean=rep(0.05, 2), sigma=diag(2) * 0.05)

hotelling <- cset(dat=bivar, method="hotelling", alpha=0.1)
summary(hotelling)
plot(hotelling, main="90\% Hotelling Region")

# Example 2: simultaneous 90\% confidence region for the mean and variance of univariate normal data

univar <- rnorm(n=50)

moodvar <- csetMV(dat=univar, method="mood", alpha=0.1, scale="var")
summary(moodvar)
plot(moodvar, main="90\% Mood Region")
}
}

\keyword{classes}