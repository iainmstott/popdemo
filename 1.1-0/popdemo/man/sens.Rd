% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sens.R
\name{sens}
\alias{sens}
\title{Calculate sensitivity matrix}
\usage{
sens(A, eval = "max", all = FALSE)
}
\arguments{
\item{A}{a square, non-negative numeric matrix of any dimension}

\item{eval}{the eigenvalue to evaluate. Default is \code{eval="max"}, which 
evaluates the dominant eigenvalue (the eigenvalue with largest REAL value: 
for imprimitive or reducible matrices this may not be the first eigenvalue). 
Otherwise, specifying e.g. \code{eval=2} will evaluate sensitivity of the 
eigenvalue with second-largest modulus.}

\item{all}{(optional) if \code{FALSE}, then only sensitivity values for 
observed transitions (nonzero entries in \code{A}) are returned.}
}
\value{
A numeric (real or complex) matrix of equal dimension to \code{A}.
}
\description{
Calculate the sensitivity matrix for a population matrix projection model 
using eigenvectors.
}
\details{
\code{sens} uses the eigenvectors of \code{A} to calculate the sensitivity 
matrix of the specified eigenvalue, see section 9.1 in Caswell (2001). 
Same method as \code{sensitivity} in \code{popbio} but can also evaluate
subdominant eigenvalues.
}
\examples{
  # Create a 3x3 PPM
  ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )

  # Calculate sensitivities of dominant eigenvalue
  sens(A)
  # Calculate sensitivities of first subdominant eigenvalue,
  # only for observed transitions
  sens(A, eval=2, all=FALSE)

}
\references{
Caswell (2001) Matrix Population Models 2nd ed. Sinauer.
}
\seealso{
Other PerturbationAnalyses: \code{\link{elas}},
  \code{\link{tfa_inertia}}, \code{\link{tfa_lambda}},
  \code{\link{tfam_inertia}}, \code{\link{tfam_lambda}},
  \code{\link{tfs_inertia}}, \code{\link{tfs_lambda}}
}
\concept{
perturbation sensitivity elasticity
}
