\name{recovery}
\alias{recovery}
\non_function{}
\title{ Recovery Time Data Set }
\usage{data(recovery)}
\description{
 Recovery time data set taken from Westfall 
 et al. (1999, p. 80). Many-to-one comparisons according to 
 Dunnett in an unbalanced one-way layout. 
}
\format{
  This data frame contains the following variables
  \describe{
    \item{blanket}{Blanket type at 4 levels: \code{b0}, \code{b1},
                   \code{b3}, \code{b4}.}
    \item{minutes}{Response variable: recovery time after a surgical 
                   procedure.}
  }
}
\details{
  See Westfall et al. (1999, p. 80)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{

data(recovery)

amod <- aov(minutes ~ blanket, data = recovery)

### cf. Westfall et al. (1999, p. 80)
confint(mcp(amod, hypotheses = list(blanket = "Dunnett"), alternative = "less"), level = 0.9)

}
\keyword{datasets}
