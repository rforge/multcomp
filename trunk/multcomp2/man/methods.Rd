\name{glht-methods}
\alias{summary.glht}
\alias{confint.glht}
\alias{coef.glht}
\alias{vcov.glht}
\alias{univariate}
\alias{adjusted}
\alias{global}
\title{ General Linear Hypotheses }
\description{
  Simultaneous Inference for General Linear Hypotheses
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
  \item{object}{an object of class \code{\link{glht}}.}
  \item{test}{a function for computing p values.}
  \item{parm}{}
  \item{level}{confidence level}
  \item{null}{logical, indicating whether 
              \eqn{K \hat{\beta}} or \eqn{m} (\code{null = FALSE})
              should be returned.}
  \item{type}{the type of adjustment to be applied.}
  \item{...}{additional arguments, such as \code{maxpts}, 
             \code{abseps} or \code{releps} to
             \code{\link[mvtnorm]{pmvnorm}} in \code{adjusted} or 
             \code{\link[mvtnorm]{qmvnorm}} in \code{confint}.}
}
\details{

  

}
\value{

  \code{summary} computes (adjusted) p values for general linear hypotheses,
  \code{confint} computes adjusted confidence intervals. 
  \code{coef} returns the estimates linear combinations and \code{vcov}
  its covariance. 

}
\examples{

  ### set up a two-way ANOVA with interactions
  amod <- aov(breaks ~ wool * tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  wht <- glht(amod, linfct = mcp(tension = "Tukey"))

  ### 95% simultaneous confidence intervals
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
