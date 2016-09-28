\name{ Cohen.cumulative }
\alias{ Cohen.cumulative }
\title{ Calculate Cohen's cumulative distance }
\description{
    Calculate Cohen's cumulative distance metric for a
    given population projection matrix (PPM) model.
}
\usage{
    Cohen.cumulative(A, vector)
}
\arguments{
  \item{A}{
    a square, irreducible, non-negative numeric matrix of any dimension.
  }
  \item{vector}{
    a specified initial age/stage distribution of class vector or 
    class matrix.
  }
}
\details{
    Calculates the cumulative distance metric as outlined in Cohen (1979).\cr

    \code{Cohen.cumulative} will not work for reducible matrices and returns a warning
    for imprimitive matrices (although will not function for imprimitive matrices with
    nonzero imaginary components in the dominant eigenpair).
}
\value{
    Cohen's D1.
}
\references{
    Cohen (1979) SIAM J. Appl. Math., 36, 169-175.\cr
    Stott et al. (2011) Ecol. Lett., 14, 959-970.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate Cohen cumulative distance
    Cohen.cumulative(A, vector=initial)
}
\keyword{ Cohen }
\keyword{ cumulative }
\keyword{ distance }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }