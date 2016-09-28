\name{ is.matrix_ergodic }
\alias{ is.matrix_ergodic }
\title{ Determine ergodicity of a matrix }
\description{
    Determine whether a matrix is ergodic
    or nonergodic
}
\usage{
    is.matrix_ergodic(A, digits=12, return.eigvec=FALSE)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension.
  }
  \item{digits}{
    the number of digits that the dominant left eigenvector should be rounded to.
  }
  \item{return.eigvec}{
    (optional) logical argument determining whether or not the dominant left eigenvector
    should be returned.
  }
}
\details{
    \code{is.matrix_ergodic} works on the premise that a matrix is ergodic if and only if
    the dominant left eigenvector (the reproductive value vector) of the matrix is positive 
    (Stott et al. 2010).\cr

    In rare cases, \R may calculate that the dominant left eigenvector of a nonergodic
    matrix contains very small entries that are approximate to (but not equal to) zero.  Rounding
    the dominant eigenvector using \code{digits} prevents mistakes.
}
\value{
    If \code{return.eigvec=FALSE}, either \code{TRUE} (for an ergodic matrix) or \code{FALSE}
    (for a nonergodic matrix).\cr
    If \code{return.eigvec=TRUE}, a list containing elements:\cr
  \item{ergodic}{
    \code{TRUE} or \code{FALSE}
  }
  \item{eigvec}{
    the dominant left eigenvector of \code{A}.
  }
}
\references{
    Stott et al. (2010) Methods Ecol. Evol., 1, 242-252.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{is.matrix_primitive}}, \code{\link{is.matrix_irreducible}}.
}
\examples{
    # Create a 3x3 ergodic PPM
    A <- matrix(c(0,0,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Diagnose ergodicity
    is.matrix_ergodic(A)

    # Create a 3x3 nonergodic PPM
    B<-A
    B[3,2] <- 0
    B

    # Diagnose ergodicity and return
    # left eigenvector
    is.matrix_ergodic(B, return.eigvec=TRUE)
}
\keyword{ ergodicity }
\keyword{ ergodic }
\keyword{ nonergodic }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }