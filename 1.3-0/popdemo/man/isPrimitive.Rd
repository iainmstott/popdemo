% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isPrimitive.R
\name{isPrimitive}
\alias{isPrimitive}
\title{Determine primitivity of a matrix}
\usage{
isPrimitive(A)
}
\arguments{
\item{A}{a square, non-negative numeric matrix of any dimension.}
}
\value{
\code{TRUE} (for an primitive matrix) or \code{FALSE} (for an imprimitive 
matrix).
}
\description{
Determine whether a matrix is primitive or imprimitive
}
\details{
\code{isPrimitive} works on the premise that a matrix \strong{A} is 
primitive if \strong{A}^(s^2-(2*s)+2) is positive, where s is the dimension 
of \strong{A} (Caswell 2001).
}
\examples{
  # Create a 3x3 primitive PPM
  ( A <- matrix(c(0,1,2,0.5,0,0,0,0.6,0), byrow=TRUE, ncol=3) )

  # Diagnose primitivity
  isPrimitive(A)

  # Create a 3x3 imprimitive PPM
  B<-A; B[1,2] <- 0; B

  # Diagnose primitivity
  isPrimitive(B)

}
\references{
\itemize{
 \item Caswell (2001) matrix Population Models, 2nd. ed. Sinauer.
}
}
\seealso{
Other PerronFrobeniusDiagnostics: \code{\link{isErgodic}},
  \code{\link{isIrreducible}}
}
\concept{
primitivity primitive imprimitive Perron Frobenius
}
