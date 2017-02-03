\name{marzo}

\alias{marzo}

\docType{data}

\title{
Bioequivalence study of ticlopidine hydrochloride
}

\description{
Pharmacokinetic data from a study on the bioequivalence of a test and a reference formulation of ticlopidine hydrochloride in 24 healthy male volunteers, using a randomised crossover design (Marzo et al. 2002).
}

\usage{data("marzo")}

\format{
  A data frame with 24 observations on the following 8 variables.
  \describe{
    \item{\code{Volunteer}}{A numeric vector giving the volunteer ID.}
    \item{\code{Sequence}}{A factor with levels \code{RT} and \code{TR} specifying the sequence a volunteer was randomised to (R=reference, T=test).}
    \item{\code{Cmax_T}}{A numeric vector of the maximum concentration (Cmax) with the test product.}
    \item{\code{Cmax_R}}{A numeric vector of the maximum concentration (Cmax) with the reference product.}
    \item{\code{AUC_T}}{A numeric vector of the area under the concentration-time curve (AUC) from zero to the last observed time point with the test product.}
    \item{\code{AUC_R}}{A numeric vector of the area under the concentration-time curve (AUC) from zero to the last observed time point with the reference product.}
    \item{\code{AUCinf_T}}{A numeric vector of the area under the concentration-time curve AUC) from zero to infinity with the test product.}
    \item{\code{AUCinf_R}}{A numeric vector of the area under the concentration-time curve AUC) from zero to infinity with the reference product.}
  }
}

\details{
The pharmacokinetic parameters (Cmax and AUC) were calculated using a non-compartmental approach. The data were taken from Tables I and II of Marzo et al. (2002).
}

\source{
Antonio Marzo, Lorenzo Dal Bo, Antonio Rusca, Pierangelo Zini (2002) Bioequivalence of ticlopidine hydrochloride administered in single dose to healthy volunteers. Pharmacological Research, 46(5), 401--407.
}

\references{
Philip Pallmann & Thomas Jaki (2017) Simultaneous confidence regions and intervals for multivariate bioequivalence. Submitted to Statistics in Medicine.
}

\examples{
data(marzo)

## An example analysis of Cmax assuming log-normality

# Difference of log(Cmax)
marzo$deltalogCmax <- log(marzo$Cmax_T) - log(marzo$Cmax_R)

# Estimated mean treatment effect with SE 
mean(marzo$deltalogCmax)
sd(marzo$deltalogCmax) / sqrt(nrow(marzo))

# Two one-sided test (TOST) p-values
t.test(x=marzo$deltalogCmax, alternative="less", mu=log(1.25))$p.value
t.test(x=marzo$deltalogCmax, alternative="greater", mu=log(0.80))$p.value

# 90\% TOST confidence interval
t.test(x=marzo$deltalogCmax, conf.level=0.9)$conf.int[1:2]
}

\keyword{datasets}
