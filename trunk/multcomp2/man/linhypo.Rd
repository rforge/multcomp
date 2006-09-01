\name{linhypo}
\alias{linhypo}
\title{ Linear Hypotheses }
\description{
  Compute Matrices of Linear Hypotheses for Multiple Comparison Procedures 
}
\usage{
linhypo(x, type = c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", 
                    "Williams", "Marcus", "McDermott"), base = 1)
}
\arguments{
  \item{x}{ a factor giving the groups to be compared. }
  \item{type}{ a character string specifying the type to linear hypotheses to be computed
               for the levels of \code{x}.}
  \item{base}{ an integer specifying which group is considered the baseline
          group for \code{type = "Dunnett"}.}
}
\details{

  

}
\value{
  A matrix of linear hypotheses. 
}
\examples{

  ### Tukey all-pair comparisons
  linhypo(gl(3, 9), "Tukey")

}
\keyword{misc}
