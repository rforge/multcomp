\name{glht}
\alias{glht}
\title{ General Linear Hypotheses }
\description{
  General Linear Hypotheses in the Generalized Linear Model
}
\usage{
glht(model, K, m = 0, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{model}{a fitted (generalized) linear model, 
               for example \code{\link{lm}}, 
               \code{\link{glm}}, \code{\link{aov}} etc.}
  \item{K}{ either a matrix of linear hypotheses for 
            the coefficients of \code{model} or a symbolic
            description of the linear hypotheses. Alternatively,
            a named list of matrices, symbolic descriptions 
            or \code{type} arguments for \code{\link{contrMat}} 
            for multiple comparisons of means, i.e., 
            for the levels of the corresponding factors. }
  \item{m}{right hand side of the hypothesis: \eqn{H_0: K \beta = m}.}
  \item{alternative}{a character string specifying the alternative 
      hypothesis, must be one of '"two.sided"' (default), '"greater"' or
      '"less"'.  You can specify just the initial letter.}
}
\details{

  

}
\value{

  An object of class \code{glht} with \code{summary}, 
  \code{confint}, \code{coef} and 
  \code{vcov} methods being available.

}
\examples{

  ### set up a two-way ANOVA with interactions
  amod <- aov(breaks ~ wool * tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  glht(amod, K = list(tension = "Tukey"))


  ### simultaneous p values for parameters in a Cox model
  leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), data = leuk)

  ### set up linear hypotheses
  glht(leuk.cox, K = diag(length(coef(leuk.cox))))

}
\keyword{htest}
