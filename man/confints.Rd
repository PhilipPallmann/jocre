\name{confints}

\alias{confints}

\title{Confidence intervals around multivariate normal means}

\description{Computes boundaries of (simultaneous) confidence intervals around multivariate normal means using different methods.}

\usage{confints(dat, method, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{A matrix or data.frame with independent units in rows and multivariate outcomes in columns.}
  \item{method}{A character string specifying the method to be used. See details for available \code{method}s.}
  \item{alpha}{A numeric value giving the type I error level to be controlled. Default is \code{0.1}.}
  \item{steps}{An integer setting the initial number of steps for the search algorithm.}
}

\details{Available \code{method}s are: \code{expanded} for the two one-sided test (TOST) procedure (Schuirmann 1987) using the expanded intervals described e.g., in Bofinger (1992) and Hsu et al. (1994); \code{fixseq} for the fixed sequence intervals described in Maurer et al (1995) and Hsu & Berger (1999); \code{tost} for the two one-sided test (TOST) intervals described in Schuirmann (1987).}

\value{A matrix of parameter estimates and confidence bounds.}

\references{
Eve Bofinger (1992) Expanded confidence intervals, one-sided tests, and equivalence testing. Journal of Biopharmaceutical Statistics, 2(2), 181--188.

Jason C. Hsu & Roger L. Berger (1999) Stepwise confidence intervals without multiplicity adjustment for dose-response and toxicity studies. Journal of the American Statistical Association, 94(446), 468--482.

Jason C. Hsu, J. T. Gene Hwang, Hung-Kung Liu, Stephen J. Ruberg (1994) Confidence intervals associated with tests for bioequivalence. Biometrika, 81(1), 103--114.

Willi Maurer, Ludwig A. Hothorn, Walter Lehmacher (1995) Multiple comparisons in drug clinical trials and preclinical assays: a priori ordered hypotheses. In: Joachim Vollmar (editor), Biometrie in der Chemisch-Pharmazeutischen Industrie, vol. 6, pp. 3--18. Fischer-Verlag, Stuttgart, Germany.

Philip Pallmann & Thomas Jaki (2016) Simultaneous confidence regions and intervals for multivariate bioequivalence. Submitted to Statistics in Medicine.

Donald J. Schuirmann (1987) A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and Biopharmaceutics, 15(6), 657--680.
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

\note{Warning: please use with care! Some of the functionality has not yet been thoroughly tested.}

\seealso{\code{\link{plot2D}} for plots.}

\examples{
# Generate trivariate normal data

mydata <- mvtnorm::rmvnorm(n=50, mean=rep(0.1, 3), sigma=diag(3) * 0.05)
colnames(mydata) <- LETTERS[1:3]

# Compute 90% confidence intervals

confints(dat=mydata, method="tost", alpha=0.1)
confints(dat=mydata, method="expanded", alpha=0.1)
#confints(dat=mydata, method="fixseq", alpha=0.1)
}

\keyword{multivariate}
