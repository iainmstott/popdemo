\name{ projection.distance }
\alias{ projection.distance }
\title{ Calculate projection distance }
\description{
    Calculate projection distance for a given population
    projection matrix (PPM) model.
}
\usage{
    projection.distance(A, vector)
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
    \code{projection.distance} (Haridas & Tuljapurkar 2007) is the difference between
    the reproductive value of a population with demographic distribution given
    by \code{vector} and the reproductive value of a population in stable state.\cr

    \code{projection.distance} will not work for reducible matrices and returns a warning
    for imprimitive matrices (although will not function for imprimitive matrices with
    nonzero imaginary components in the dominant eigenpair).
}
\value{
    Projection distance.
}
\references{
    Haridas & Tuljapurkar (2007) Ecol. Lett., 10, 1143-1153.\cr    
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

    # Calculate projection distance
    projection.distance(A, vector=initial)
}
\keyword{ projection distance }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }