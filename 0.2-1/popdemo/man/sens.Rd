\name{ sens }
\alias{ sens }
\title{ Calculate sensitivity matrix }
\description{
    Calculate the sensitivity matrix using eigenvectors
    for a specified population projection matrix (PPM).
}
\usage{
    sens(A, eval="max", all=FALSE)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{eval}{
    the eigenvalue to evaluate.  Default is \code{eval="max"}, which evaluates the dominant
    eigenvalue (the eigenvalue with largest REAL value: for a reducible matrix this may not be the 
    first eigenvalue).  Otherwise, specifying e.g. \code{eval=2} will evaluate sensitivity of the
    eigenvalue with second-largest modulus.
  }
  \item{all}{
    (optional) if \code{FALSE}, then only sensitivity values for observed transitions are returned.
  }
}
\details{
    \code{sens} uses the eigenvectors of \code{A} to calculate the sensitivity matrix of the
    specified eigenvalue, see section 9.1 in Caswell (2001).
}
\value{
    A sensitivity matrix of equal dimension to \code{A}.
}
\references{
    Caswell (2001) Matrix Population Models 2nd ed. Sinauer.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{elas}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)

    # Calculate sensitivities of dominant eigenvalue
    sens(A)

    # Calculate sensitivities of first subdominant eigenvalue,
    # only for observed transitions
    sens(A, eval=2, all=FALSE)
}
\keyword{ sensitivity }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }