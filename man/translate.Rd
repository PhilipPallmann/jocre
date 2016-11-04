\name{translate}

\alias{translate}

\title{Conversion of confidence levels between intervals and ellipses}

\description{Translates the confidence level of a joint 100(1 -- alpha)\% confidence ellipse into that of the corresponding marginal confidence interval when projecting the ellipse's boundary onto the axes. Also does the "inverse operation" i.e., calculates the confidence level of a joint confidence ellipse so that its perpendicular shadows onto the axes are 100(1 -- alpha)\% confidence intervals.}

\usage{
translate(level=0.95, ddf, direction)
}

\arguments{
  \item{level}{A numeric value giving the confidence level.}
  \item{ddf}{An integer specifying the denominator degrees of freedom. Setting this to \code{0} enforces an asymptotic computation.}
  \item{direction}{A character string indicating what is to be computed. Choose either \code{ci2cr} or \code{cr2ci} (see details).}
}

\details{
Setting \code{direction="ci2cr"} calculates the confidence level of a confidence interval generating ellipse (CIGE) whose perpendicular shadows onto the axes are 100(1 -- alpha)\% confidence intervals with a marginal confidence level (1 -- alpha) as specified in \code{level}; see p. 205 of Fox (2008).

On the other hand, setting \code{direction="cr2ci"} computes the marginal confidence level of the intervals obtained by projecting a joint 100(1 -- alpha)\% confidence ellipse with (1 -- alpha) as specified in \code{level}; see p. 254 of Monette (1990). These marginal intervals can be viewed as including a Scheffe penalty (Scheffe 1953).

For \code{ddf=0} the F-distribution used for calculating the confidence levels is replaced with an asymptotic chi-square distribution.
}

\value{A numeric value giving the calculated confidence level.}

\references{
John Fox (2008) Applied Linear Regression and Generalized Linear Models. Second Edition. SAGE, Thousand Oaks, CA.

Georges Monette (1990). Geometry of multiple regression and interactive 3-D graphics. In: John Fox & J. Scott Long (eds.) Modern Methods of Data Analysis. SAGE, Newbury Park, CA.

Henry Scheffe (1953) A method for judging all contrasts in the analysis of variance. Biometrika, 40(1--2), 87--104.
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

\examples{
# Get CIGE level
translate(0.95, ddf=1, "ci2cr")
translate(0.95, ddf=9999, "ci2cr")
translate(0.95, ddf=0, "ci2cr")

# Get Scheffe CI level
translate(0.95, ddf=1, "cr2ci")
translate(0.95, ddf=9999, "cr2ci")
translate(0.95, ddf=0, "cr2ci")
}

\keyword{multivariate}
