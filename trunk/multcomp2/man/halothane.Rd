\name{halothane}
\alias{halothane}
\non_function{}
\title{ Repeated Measures Data Set }
\usage{data(Halothane)}
\description{
 Repeated measures on dogs data set taken from Westfall et al. (1999, p. 214). 2x2 factorial application of Halothane and CO2 pressure as repeated measurements on dogs with unstructured covariance matrix between the subjects and Satterthwaite correction for the df.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{treatment}{a factor at four levels: \code{HA}, \code{LA}, \code{HP}, and \code{LP}.}
    \item{dog}{Subject ID at 19 levels: 1 ... 19}
    \item{rate}{Response variable: milliseconds between heartbeats}
  }
}
\details{
  See Westfall et al. (1999, p. 214)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{

  data(halothane)
  amod <- aov(rate ~ treatment + dog, data=halothane)
  confint(mcp(amod, hypotheses = list(treatment = "Tukey")))

}
\keyword{datasets}
