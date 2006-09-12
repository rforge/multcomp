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
               for example an object returned by \code{\link{lm}}, 
               \code{\link{glm}}, \code{\link{aov}} etc. \code{\link{coef}},
               \code{\link{vcov}}, \code{\link{model.matrix}}, 
               \code{\link{model.frame}} and \code{\link{terms}}
               methods are expected to be available for \code{model}.}
  \item{K}{ a specification of the linear hypotheses to be tested. 
            \code{K} can be supplied in form of a matrix of coefficients or
            by a symbolic description interpretable as expression.
            Alternatively, \code{K} may be a named list of matrices, 
            symbolic descriptions or \code{type} arguments 
            for \code{\link{contrMat}} when multiple comparisons of 
            means, i.e., factor levels, are of special interest.}
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
