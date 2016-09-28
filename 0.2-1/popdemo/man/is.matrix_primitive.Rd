\name{ is.matrix_primitive }
\alias{ is.matrix_primitive }
\title{ Determine primitivity of a matrix }
\description{
    Determine whether a matrix is primitive
    or imprimitive
}
\usage{
    is.matrix_primitive(A)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension.
  }
}
\details{
    \code{is.matrix_primitive} works on the premise that a matrix A is primitive if A raised to the power 
    of s^2-(2*s)+2 is positive, where s is the dimension of A (see Caswell 2001).
}
\value{
    \code{TRUE} (for an primitive matrix) or \code{FALSE} (for an imprimitive matrix).
}
\references{
    Caswell (2001) matrix Population Models, 2nd. ed. Sinauer.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{is.matrix_irreducible}}, \code{\link{is.matrix_ergodic}}.
}
\examples{
    # Create a 3x3 primitive PPM
    A <- matrix(c(0,1,2,0.5,0,0,0,0.6,0), byrow=TRUE, ncol=3)
    A

    # Diagnose primitivity
    is.matrix_primitive(A)

    # Create a 3x3 imprimitive PPM
    B<-A
    B[1,2] <- 0
    B

    # Diagnose primitivity
    is.matrix_primitive(B)
}
\keyword{ primitivity }
\keyword{ primitive }
\keyword{ imprimitive }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }