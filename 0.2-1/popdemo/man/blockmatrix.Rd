\name{ blockmatrix }
\alias{ blockmatrix }
\title{ Block-permute a reducible matrix }
\description{
    Block-permute a reducible matrix
}
\usage{
    blockmatrix(A)
}
\arguments{
  \item{A}{
    a square, reducible, non-negative numeric matrix of any dimension
  }
}
\details{
    Any reducible matrix may have its rows and columns simultaneously permuted so that
    it takes a specific structure, with irreducible square submatrices on the diagonal, zero
    submatrices in the lower triangle and nonzero submatrices in the upper triangle (see
    Caswell 2001; Stott et al. 2010).  \code{blockmatrix} permutes a reducible matrix into this form,
    which enables further evaluation (e.g. computation of eigenvalues of submatrices).
}
\value{
    A list containing components:\cr
  \item{blockmatrix}{    
    the block-permuted matrix.
  }
  \item{order}{
    the permutation of rows/columns of \code{A} in the block-permuted matrix.
  }
}
\references{
    Caswell (2001) Matrix population models 2nd ed. Sinauer.\cr
    Stott et al. (2010) Methods. Ecol. Evol., 1, 242-252.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\examples{
    # Create a 3x3 reducible PPM
    A <- matrix(c(0,1,0,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Block-permute the matrix
    blockmatrix(A)
}
\keyword{ reducible }
\keyword{ irreducible }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }