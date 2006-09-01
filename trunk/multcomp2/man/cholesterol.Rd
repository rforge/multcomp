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
  are defined as: 20mg once a day (\code{1time}), 10mg twice a day
  (\code{2times}) and 5mg four times a day (\code{4times}). 
  In addition, two competing drugs were used as control: (\code{drugD}) 
  and (\code{drugE}).

}
\source{

  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  

}
\examples{

  data("cholesterol")

  ### adjusted p-values for all-pairwise comparisons in a one-way layout 
  ### set up ANOVA model  
  amod <- aov(response ~ trt, data = cholesterol)

  ### set up multiple comparisons object for all-pair comparisons
  cmcp <- mcp(amod, hypotheses = list(trt = "Tukey"))

  ### cf. Westfall et al. (1999, page 171)
  summary(cmcp, distribution = univariate())
  summary(cmcp, distribution = adjusted("Shaffer"))
  summary(cmcp, distribution = adjusted("Westfall"))

  ### use only a subset of all pairwise hypotheses
  K <- linhypo(cholesterol$trt, type="Tukey")
  Ksub <- rbind(K[c(1,2,5),],
                "D-test" = c(-1, -1, -1, 3, 0),
                "E-test" = c(-1, -1, -1, 0, 3))

  ### reproduce results in Westfall et al. (1999, page 172)
  amod <- aov(response ~ trt - 1, data = cholesterol)
  summary(mcp(amod, hypotheses = list(trt = Ksub[,5:1])), 
          distribution = adjusted("Westfall"))

}
\keyword{datasets}
