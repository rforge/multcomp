\name{litter}
\alias{litter}
\non_function{}
\title{ Litter Weights Data Set }
\usage{data(litter)}
\description{
  Dose response of litter weights in rats; data set taken from Westfall et al. (1999, p. 109). Specific dose comparisons in a dose response experiment with two covariates. 
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
}
\examples{
data(litter)

amod <- aov(weight ~ dose + gesttime + number, data = litter)

# define matrix of linear hypotheses for `dose'
K <- rbind(c(1, -1,  0,  0),
           c(1,  0, -1,  0),
           c(1,  0,  0, -1),
           c(0.75, 0.25, -0.25, -0.75),
           c(0.384, 0.370, 0.246, -1),
           c(0.887, 0.113, -0.339, -0.661))
rownames(K) <- c("cont-low", "cont-mid", "cond-high", "ordinal", "arith", "log ord")

### cf. Westfall et al. (1999, p. 111)
summary(mcp(amod, hypotheses = list(dose = K), alternative = "greater"))
confint(mcp(amod, hypotheses = list(dose = K), alternative = "greater"))

}
\keyword{datasets}
