\name{ is.matrix_irreducible }
\alias{ is.matrix_irreducible }
\title{ Determine reducibility of a matrix }
\description{
    Determine whether a matrix is irreducible
    or reducible
}
\usage{
    is.matrix_irreducible(A)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension.
  }
}
\details{
    \code{is.matrix_irreducible} works on the premise that a matrix is irreducible if and only if
    (I+A) raised to the power of s-1 is positive, where I is the identity matrix of the same
    dimension as \code{A} and s is the dimension of A (see Caswell 2001).
}
\value{
    \code{TRUE} (for an irreducible matrix) or \code{FALSE} (for a reducible matrix).
}
\references{
    Caswell (2001) matrix Population Models, 2nd. ed. Sinauer.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{is.matrix_primitive}}, \code{\link{is.matrix_ergodic}}.
}
\examples{
    # Create a 3x3 irreducible PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Diagnose reducibility
    is.matrix_irreducible(A)

    # Create a 3x3 reducible PPM
    B<-A
    B[3,2] <- 0
    B

    # Diagnose reducibility
    is.matrix_irreducible(B)
}
\keyword{ reducibility }
\keyword{ irreducible }
\keyword{ reducible }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }