\name{glht-methods}
\alias{summary.glht}
\alias{confint.glht}
\alias{coef.glht}
\alias{vcov.glht}
\alias{univariate}
\alias{adjusted}
\alias{global}
\title{ Methods for General Linear Hypotheses }
\description{
  Simultaneous tests and confidence intervals for general linear 
  hypotheses.
}
\usage{
\method{summary}{glht}(object, test = adjusted(), ...)
\method{confint}{glht}(object, parm, level = 0.95, ...)
\method{coef}{glht}(object, null = FALSE, ...)
\method{vcov}{glht}(object, ...)
univariate()
adjusted(type = c("free", "Bonferroni", "Shaffer", "Westfall"), ...)
global(type = c("Chisq", "F"))
}
\arguments{
  \item{object}{ an object of class \code{\link{glht}}.}
  \item{test}{ a function for computing p values.}
  \item{parm}{ additional parameters, currently ignored.}
  \item{level}{ the confidence level required.}
  \item{null}{logical, indicating whether the linear function
              \eqn{K \hat{\beta}} or the right hand side 
              \eqn{m} (\code{null = FALSE}) of the linear hypothesis
              should be returned.}
  \item{type}{ the type of test statistic (\code{global}) or 
               multiplicity adjustment (\code{adjusted}) 
               to be applied. See below.}
  \item{...}{ additional arguments, such as \code{maxpts}, 
              \code{abseps} or \code{releps} to
              \code{\link[mvtnorm]{pmvnorm}} in \code{adjusted} or 
              \code{\link[mvtnorm]{qmvnorm}} in \code{confint}.}
}
\details{

  The methods for general linear hypotheses as described by objects returned
  by the \code{\link{glht}} can be used to actually test the global
  null hypothesis, each of the partial hypotheses and for
  simultaneous confidence intervals for the linear function $K \beta$.
  
  The \code{coef} and \code{vcov} methods compute the linear
  function $K \hat{\beta}$ and its covariance, respectively.

  The \code{test} argument to \code{summary} takes a function specifying
  the type of test to be applied. Classical Chisq or F statistics
  for testing the global hypothesis $H_0$ are implemented in function
  \code{global}. Multiplicity adjusted p values for each of the 
  linear hypotheses are implemented in function \code{adjusted}.
  
  Different multiplicity adjustments are available. The \code{type}
  argument to \code{adjusted} specifies the method to be applied:
  \code{"free"} implements adjusted p values based on the joint
  normal or $t$ distribution of the linear function, \code{"Bonferroni"}
  implements a simple Bonferroni-adjustment, and \code{"Shaffer"}
  and \code{"Westfall"} implements logically constraint 
  multiplicity adjustments.

  Simultanous confidence intervals for linear functions can be computed
  using method \code{confint}.

  All simultaneous inference procedures implemented here control
  the family-wise error rate (FWER). Multivariate
  normal and $t$ distributions are evaluated using the procedures
  implemented in package \code{mvtnorm}.

}
\value{

  \code{summary} computes (adjusted) p values for general linear hypotheses,
  \code{confint} computes adjusted confidence intervals. 
  \code{coef} returns the estimates linear combinations and \code{vcov}
  its covariance. 

}
\references{

   Juliet P. Shaffer (1986), 
   Modified sequentially rejective multiple test procedures. 
   \emph{Journal of the American Statistical Association},
   \bold{81}, 826--831.

   Peter H. Westfall (1997), 
   Multiple testing of general contrasts using logical constraints
   and correlations. \emph{Journal of the American Statistical Association},
   \bold{92}, 299--306.

}
\examples{

  ### set up a two-way ANOVA with interactions
  amod <- aov(breaks ~ wool * tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  wht <- glht(amod, linfct = mcp(tension = "Tukey"))

  ### 95\% simultaneous confidence intervals
  plot(print(confint(wht)))

  ### the same (for balanced designs only)
  TukeyHSD(amod, "tension")

  ### corresponding adjusted p values
  summary(wht)


  ### confidence bands for a simple linear model, `cars' data
  plot(cars, xlab = "Speed (mph)", ylab = "Stopping distance (ft)",
       las = 1)

  ### fit linear model and add regression line to plot
  lmod <- lm(dist ~ speed, data = cars)
  abline(lmod)

  ### a grid of speeds
  speeds <- seq(from = min(cars$speed), to = max(cars$speed), 
                length = 10)

  ### linear hypotheses: 10 selected points on the regression line != 0
  K <- cbind(1, speeds)                                                        

  ### set up linear hypotheses
  cht <- glht(lmod, linfct = K)

  ### confidence intervals, i.e., confidence bands, and add them plot
  cci <- confint(cht)
  lines(speeds, cci$confint[,"lwr"], col = "blue")
  lines(speeds, cci$confint[,"upr"], col = "blue")


  ### simultaneous p values for parameters in a Cox model
  if (require("survival") && require("MASS")) {
      data("leuk", package = "MASS")
      leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), data = leuk)

      ### set up linear hypotheses
      lht <- glht(leuk.cox, linfct = diag(length(coef(leuk.cox))))

      ### adjusted p values
      print(summary(lht))
  }

}
\keyword{htest}
