\name{mmm}

\alias{mmm}
\alias{mlf}

\title{Simultaneous Inference for Multiple Marginal Models}

\description{Calculation of correlation between test statistics from multiple marginal models 
             using the score decomposition}

\usage{
mmm(...)
mlf(...)
}

\arguments{
  \item{\dots}{A names argument list containing fitted models (\code{mmm}) or
               definitions of linear functions (\code{mlf}).}
}

\details{

  Estimated correlations of the estimated parameters of interest from the
  multiple marginal models are obtained using a stacked version of the i.i.d. 
  decomposition of parameter estimates by means of score components (first
  derivatives of the log likelihood).  The method is less conservative than
  the Bonferroni correction.  The details are provided by Pipper and Ritz
  (2012).

}

\value{

  An object of class \code{mmm} or \code{mlf}, basically a named list of the
  arguments with a special method for \code{\link{glht}} being available for
  the latter.  \code{vcov}, \code{\link[sandwich]{estfun}}, and
  \code{\link[sandwich]{bread}} methods are available for objects of class
  \code{mmm}.

}

\references{
  
  Pipper, C. B., Ritz, C. and Bisgaard, H. (2012) A versatile method for confirmatory evaluation of the effects of a covariate in multiple models,
  \emph{Journal of the Royal Statistical Society, Series C Applied Statistics}, \bold{61}, 315--326.
}

\author{

  Code for the computation of the joint covariance and sandwich matrix was contributed by 
  Christian Ritz and Christian B. Pipper.

}

\examples{

## Combining linear regression and logistic regression
set.seed(29)
y1 <- rnorm(100)
y2 <- factor(y1 + rnorm(100, sd = .1) > 0)
x1 <- gl(4, 25) 
x2 <- runif(100, 0, 10)

m1 <- lm(y1 ~ x1 + x2)
m2 <- glm(y2 ~ x1 + x2, family = binomial())
## Note that the same explanatory variables are considered in both models
## but the resulting parameter estimates are on 2 different scales (original and log-odds scales)

## Simultaneous inference for the same parameter in the 2 model fits
summary(glht(mmm(m1 = m1, m2 = m2), mlf(m1 = "x12 = 0", m2 = "x12 = 0")))

## Simultaneous inference for different parameters in the 2 model fits
summary(glht(mmm(m1 = m1, m2 = m2),
             mlf(m1 = "x12 = 0", m2 = "x13 = 0")))

## Simultaneous inference for different and identical parameters in the 2 model fits
summary(glht(mmm(m1 = m1, m2 = m2),
             mlf(m1 = c("x12 = 0", "x13 = 0"), m2 = "x13 = 0")))

## Examples for binomial data

## Two independent outcomes
y1.1 <- rbinom(100, 1, 0.5)
y1.2 <- rbinom(100, 1, 0.5)
group <- factor(rep(c("A", "B"), 50))

m1 <- glm(y1.1 ~ group, family = binomial)
m2 <- glm(y1.2 ~ group, family = binomial)

summary(glht(mmm(m1 = m1, m2 = m2), mlf(m1 = "groupB = 0", m2 = "groupB = 0")))

## Two perfectly correlated outcomes
y2.1 <- rbinom(100, 1, 0.5)
y2.2 <- y2.1
group <- factor(rep(c("A", "B"), 50))

m1 <- glm(y2.1 ~ group, family = binomial)
m2 <- glm(y2.2 ~ group, family = binomial)

summary(glht(mmm(m1 = m1, m2 = m2), mlf(m1 = "groupB = 0", m2 = "groupB = 0")))

}
\keyword{models}

