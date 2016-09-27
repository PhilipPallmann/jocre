\name{iutsize}
\alias{iutsize}

\title{Size of intersection union tests derived from confidence sets}

\description{Computes the actual size of a intersection union test procedure that corresponds to a (1 -- alpha) confidence set.}

\usage{iutsize(p, n, alpha=0.1, sim=1e6)}

\arguments{
  \item{p}{An integer giving the number of dimensions.}
  \item{n}{An integer giving the sample size.}
  \item{alpha}{A numeric value specifying the type I error level to be controlled. Default is 0.1.}
  \item{sim}{An integer giving the number of simulations to be carried out. Default is 1 million.}
}

\details{A (1 -- alpha) confidence set can be used to derive a two one-sided tests (TOST) procedure (Schuirmann 1987) whereby type I error rate control is ensured at level alpha due to the intersection union principle (Berger 1982). The actual test size, however, is often substantially lower than alpha i.e., the approach is conservative. It is well known for the one-dimensional case that the TOST corresponding to a (1 -- alpha) confidence interval has size (1 -- alpha/2). This function computes the achieved test size with dimension \code{p} and \code{n} according to the formula on p. 399 of Wang et al (1999).}

\value{A numeric value giving the actual size of the test.}

\references{Roger L. Berger (1982) Multiparameter hypothesis testing and acceptance sampling. Technometrics, 24(4), 295--300.

Donald J. Schuirmann (1987) A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and Biopharmaceutics, 15(6), 657--680.

Weizhen Wang, J. T. Gene Hwang, Anirban DasGupta (1999) Statistical tests for multivariate bioequivalence. Biometrika, 86(2), 395--402.}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

%\note{xxx}

%\seealso{x}

\examples{
# For p=1 we get the well-known result that the 90% CI corresponds to the TOST at 5%:

#iutsize(p=1, n=20)

# With increasing dimension the test gets conservative:

#iutsize(p=2, n=20)
#iutsize(p=3, n=20)

# For p>1 the conservativeness also depends on sample size:

#iutsize(p=2, n=10)
#iutsize(p=2, n=1000)
}

\keyword{multivariate}
