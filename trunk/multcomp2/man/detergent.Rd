\name{detergent}
\alias{detergent}
\non_function{}
\title{ Detergent Durability Data Set }
\usage{data("detergent")}
\description{
 Detergent durability in a fractional design experiment; data set 
 taken from Westfall et al. (1999, p. 189). Balanced incomplete block design.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{detergent}{detergent at 5 levels: 1 ... 5}
    \item{block}{block at 10 levels: 1 ... 10}
    \item{plates}{Response variable: number of plates washed before 
     the foam disappears.}
  }
}
\details{
  See Westfall et al. (1999, p. 189)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{


  data("detergent")

  ### set up two-way ANOVA without interactions
  amod <- aov(plates ~ block + detergent, data = detergent)

  ### set up all-pair comparisons
  dmcp <- mcp(amod, hypotheses = list(detergent = "Tukey"))

  ### see Westfall et al. (1999, p. 190)
  confint(dmcp)

  ### see Westfall et al. (1999, p. 192)
  summary(dmcp, distribution = univariate())
  summary(dmcp, distribution = adjusted("Shaffer"))
  summary(dmcp, distribution = adjusted("Westfall"))

}
\keyword{datasets}
