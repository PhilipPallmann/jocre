\name{confset}

\alias{confset}

\title{Confidence regions around multivariate normal means}

\description{Computes boundaries of (simultaneous) confidence regions around multivariate normal means using different methods, and projects the boundaries onto the axes.}

\usage{confset(dat, method, alpha=0.1, steps=100, TsengBrownA=1, TsengBrownB=1)}

\arguments{
  \item{dat}{A matrix or data.frame with independent units in rows and multivariate outcomes in columns.}
  \item{method}{A character string specifying the method to be used. See details for available \code{method}s.}
  \item{alpha}{A numeric value giving the type I error level to be controlled. Default is 0.1.}
  \item{steps}{An integer setting the initial number of steps for the search algorithm.}
  \item{TsengBrownA}{A numeric giving the parameter A to be used when \code{method="tseng.brown"} (ignored otherwise). Default is \code{1}.}
  \item{TsengBrownB}{A numeric giving the parameter B to be used when \code{method="tseng.brown"} (ignored otherwise). Default is \code{1}.}
}

\details{Available \code{method}s are: \code{emp.bayes} for the empirical Bayes region described in Casella & Hwang (1983); \code{fixseq} for the fixed sequence intervals described in Maurer et al (1995) and Hsu & Berger (1999); \code{hotelling} for the Hotelling-type region described in Wang et al (1999); \code{limacon.asy} for the limacon-shaped mimimum expected volume region described in Brown et al (1995); \code{limacon.fin} for the finite-sample variant of the minimum expected volume region described in Berger & Hsu (1996); \code{standard.cor} for the standard region incorporating correlation between parameters; \code{standard.ind} for the standard region ignoring correlation between parameters; \code{tost} for the two one-sided test (TOST) intervals described in Schuirmann (1987); \code{tseng} for the mimimum expected interval length region described in Tseng (2002); \code{tseng.brown} for the pseudo-empirical Bayes region described in Tseng & Brown (1997). See also an overview and comparison of all methods in Pallmann & Jaki (2016).}

\value{A matrix of parameter estimates and confidence bounds.}

\references{
Roger L. Berger & Jason C. Hsu (1996) Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science, 11(4), 283--319.

Lawrence D. Brown, George Casella, J. T. Gene Hwang (1995) Optimal confidence sets, bioequivalence, and the limacon of Pascal. Journal of the American Statistical Association, 90(431), 880--889.

George Casella & Jiunn T. Hwang (1983) Empirical Bayes confidence sets for the mean of a multivariate normal distribution. Journal of the American Statistical Association, 78(383), 688--698.

%Harold Hotelling (1931) The generalization of Student's ratio. The Annals of Mathematical Statistics, 2(3), 360--378.

Jason C. Hsu & Roger L. Berger (1999) Stepwise confidence intervals without multiplicity adjustment for dose-response and toxicity studies. Journal of the American Statistical Association, 94(446), 468--482.

Willi Maurer, Ludwig A. Hothorn, Walter Lehmacher (1995) Multiple comparisons in drug clinical trials and preclinical assays: a priori ordered hypotheses. In: Joachim Vollmar (editor), Biometrie in der Chemisch-Pharmazeutischen Industrie, vol. 6, pp. 3--18. Fischer-Verlag, Stuttgart, Germany.

Philip Pallmann & Thomas Jaki (2016) Simultaneous confidence regions and intervals for multivariate bioequivalence. Submitted to Statistics in Medicine.

Donald J. Schuirmann (1987) A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and Biopharmaceutics, 15(6), 657--680.

Yu-Ling Tseng (2002) Optimal confidence sets for testing average bioequivalence. Test, 11(1), 127--141.

Yu-Ling Tseng & Lawrence D. Brown (1997) Good exact confidence sets for a multivariate normal mean. The Annals of Statistics, 25(5), 2228--2258.

Weizhen Wang, J. T. Gene Hwang, Anirban DasGupta (1999) Statistical tests for multivariate bioequivalence. Biometrika, 86(2), 395--402.

# bootkern not included so far
}

\author{Philip Pallmann (\email{p.pallmann@lancaster.ac.uk})}

\note{Warning: please use with care! Some of the functionality has not yet been thoroughly tested.}

\seealso{\code{\link{plot2D}} for plots.}

\examples{
# Generate trivariate normal data

mydata <- mvtnorm::rmvnorm(n=50, mean=rep(0.1, 3), sigma=diag(3) * 0.05)
colnames(mydata) <- LETTERS[1:3]

# Compute boundaries of simultaneous 90% confidence regions

#confset(dat=mydata, method="standard.ind", alpha=0.1)
#confset(dat=mydata, method="standard.cor", alpha=0.1)
#confset(dat=mydata, method="hotelling", alpha=0.1)
#confset(dat=mydata, method="limacon.asy", alpha=0.1)
#confset(dat=mydata, method="limacon.fin", alpha=0.1)
#confset(dat=mydata, method="tseng", alpha=0.1)

}

\keyword{multivariate}