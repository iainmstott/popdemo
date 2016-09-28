\name{ elas }
\alias{ elas }
\title{ Calculate elasticity matrix }
\description{
    Calculate the elasticity matrix using eigenanalysis
    for a specified population projection matrix (PPM).
}
\usage{
    elas(A, eval="max")
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{eval}{
    the eigenvalue to evaluate.  Default is \code{eval="max"}, which evaluates the dominant
    eigenvalue (the eigenvalue with largest REAL value: for a reducible matrix this may not be the 
    first eigenvalue).  Otherwise, specifying e.g. \code{eval=2} will evaluate elasticity of the
    eigenvalue with second-largest modulus.
  }
}
\details{
    \code{elas} uses the eigenvectors of \code{A} to calculate the elasticity matrix 
    of the specified eigenvalue, see section 9.1 in Caswell (2001).
}
\value{
    An elasticity matrix of equal dimension to \code{A}.
}
\references{
    Caswell (2001) Matrix Population Models 2nd ed. Sinauer.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{sens}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate elasticities of dominant eigenvalue
    elas(A)

    # Calculate elasticities of first subdominant eigenvalue
    elas(A, eval=2)
}
\keyword{ elasticity }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }