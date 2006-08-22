\name{cholesterol}
\alias{cholesterol}
\non_function{}
\title{ Cholesterol Reduction Data Set }
\usage{data(cholesterol)}
\description{
 Cholesterol reduction for five treatments; data set taken from 
 Westfall et al. (1999, p. 153). All pairwise comparisons according 
 to Tukey in a balanced one-way layout.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{trt}{Treatment at 5 levels, see below.}
    \item{response}{Response variable.}
  }
}
\details{
  See Westfall et al. (1999, p. 153). The treatment formulations 
are defined as: 20mg one a day (\code{1time}), 10mg twice a day
(\code{2times})
and 5mg four times a day (\code{4times}). In addition, two competing drugs
were used as control (\code{drugD}) and (\code{drugE}).
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(cholesterol)

# adjusted p-values for all-pairwise comparisons in a one-way layout 
amod <- aov(response ~ trt, data = cholesterol)
summary(mcp(amod, hypotheses = list(trt = "Tukey")))

# use only a subset of all pairwise hypotheses:
#
# * drug D versus all other formulations and
# * all pairwise comparisions for "1time", "2times" and "4times"
# the contrast matrix
K <- linhypo(cholesterol$trt, type="Tukey")
Ksub <- K[c(1,3,5,6,8,10),]
summary(mcp(amod, hypotheses = list(trt = Ksub)))

}
\keyword{datasets}
