\name{wires}

\alias{wires}

\docType{data}

\title{
Quality control study of wire connections
}

\description{
Data from a study in quality control assessing the breaking strengths of 20 wire connections between a semiconductor wafer and a terminal post (King 1971).
}

\usage{data("wires")}

\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{Strength}}{A numeric vector giving the strength at which the connection failed.}
    \item{\code{Failure}}{A factor with levels \code{b} and \code{w} specifying whether the bond or wire failed.}
  }
}

\details{
The data were taken from Table 4.1 of Nelson (1982).
}

\source{
James R. King (1971) Probability Charts for Decision Making. Industrial Press, New York, NY.
}

\references{
Wayne B. Nelson (1982) Applied Life Data Analysis. Wiley, Hoboken, NJ.
}

\examples{
data(wires)

# Simultaneous 90\% confidence regions for the mean and variance
plot(csetMV(wires$Strength, method="mood"), main="Mood")
plot(csetMV(wires$Strength, method="lrt"), main="LRT")
}

\keyword{datasets}
