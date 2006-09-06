\name{litter}
\alias{litter}
\non_function{}
\title{ Litter Weights Data Set }
\usage{data(litter)}
\description{
  Dose response of litter weights in rats; data set taken from Westfall et al. (1999, p. 109). 
  Specific dose comparisons in a dose response experiment with two covariates. 
}
\format{
  This data frame contains the following variables
  \describe{
    \item{dose}{Dosages at 4 levels: \code{0}, \code{5}, \code{50}, 
                \code{500}.}
    \item{gesttime}{Gestation time as covariate.}
    \item{number}{Number of animals in litter as covariate.}
    \item{weight}{Response variable: average post-birth weights 
                  in the entire litter.}
  }
}
\details{
  See Westfall et al. (1999, p. 109)
}
\source{

  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  

  P. H. Westfall (1997). Multiple Testing of General Contrasts Using
  Logical Constraints and Correlations. \emph{Journal of the American
  Statistical Association}, \bold{92}(437), 299--306.

}
\examples{

  data("litter")

  ### fit ANCOVA model to data
  amod <- aov(weight ~ dose + gesttime + number, data = litter)

  ### define matrix of linear hypotheses for `dose'
  doselev <- as.integer(levels(litter$dose))
  K <- rbind(contrMat(table(litter$dose), "Tukey"),
             otrend = c(-1.5, -0.5, 0.5, 1.5),
             atrend = doselev - mean(doselev),
             ltrend = log(1:4) - mean(log(1:4)))

  ### set up multiple comparison object
  Kht <- glht(amod, K = list(dose = K), alternative = "less")

  ### cf. Westfall (1997, Table 2)
  summary(Kht, distribution = univariate())
  summary(Kht, distribution = adjusted("Bonferroni"))
  summary(Kht, distribution = adjusted("Shaffer"))
  summary(Kht, distribution = adjusted("Westfall"))
  summary(Kht, distribution = adjusted("free"))

}
\keyword{datasets}
