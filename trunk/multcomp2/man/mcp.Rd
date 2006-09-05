\name{mcp}
\alias{mcp}
\alias{summary.mcp}
\alias{confint.mcp}
\alias{univariate}
\alias{adjusted}
\title{ Multiple Comparison Procedures }
\description{
  Multiple Comparison Procedures for Linear Hypotheses in the Generalized Linear Model
}
\usage{
mcp(object, hypotheses = NULL, m = 0, alternative = c("two.sided", "less", "greater"))
\method{summary}{mcp}(object, distribution = adjusted(), ...)
univariate()
adjusted(type = c("free", "Bonferroni", "Shaffer", "Westfall"), ...)
}
\arguments{
  \item{object}{ a fitted (generalized) linear model, for example \code{\link{lm}}, 
                 \code{\link{glm}}, \code{\link{aov}} etc.}
  \item{hypotheses}{ either a matrix of linear hypotheses for the coefficients of \code{object} or
                     a named list of matrices (or \code{type} arguments for \code{\link{linhypo}}) 
                     for linear hypotheses for levels of the corresponding factors. }
  \item{m}{ right hand side of the hypothesis: \eqn{H_0: K \beta = m}. }
  \item{alternative}{ a character string specifying the alternative hypothesis,
          must be one of '"two.sided"' (default), '"greater"' or
          '"less"'.  You can specify just the initial letter.}
  \item{distribution}{a function for computing p values.}
  \item{type}{the type of adjustment to be applied.}
  \item{...}{additional arguments, such as \code{maxpts}, \code{abseps} or \code{releps} to
             \code{\link[mvtnorm]{pmvnorm}} in \code{adjusted} or 
             \code{\link[mvtnorm]{qmvnorm}} in \code{confint}.}
}
\details{

  

}
\value{

  An object of class \code{mcp} with \code{summary} and \code{confint} 
  methods being available.

}
\examples{

  ### set up a two-way ANOVA with interactions
  amod <- aov(breaks ~ wool * tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  wmcp <- mcp(amod, hypotheses = list(tension = "Tukey"))

  ### 95% simultaneous confidence intervals
  plot(print(confint(wmcp)))

  ### the same (for balanced designs only)
  TukeyHSD(amod, "tension")

  ### corresponding adjusted p values
  summary(wmcp)


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

  ### set up multiple comparisons
  cmcp <- mcp(lmod, hypotheses = K)

  ### confidence intervals, i.e., confidence bands, and add them plot
  cci <- confint(cmcp)
  lines(speeds, cci$confint[,"lwr"], col = "blue")
  lines(speeds, cci$confint[,"upr"], col = "blue")


  ### simultaneous p values for parameters in a Cox model
  leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), data = leuk)

  ### set up multiple comparisons for all coefficients
  lmcp <- mcp(leuk.cox, hypotheses = diag(length(coef(leuk.cox))))

  ### adjusted p values
  summary(lmcp)

}
\keyword{htest}
