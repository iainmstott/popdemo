\name{ Keyfitz.delta }
\alias{ Keyfitz.delta }
\title{ Calculate Keyfitz's delta }
\description{
    Calculate Keyfitz's delta for a given population
    projection matrix (PPM) model.
}
\usage{
    Keyfitz.delta(A, vector)
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
    Keyfitz's delta is the sum of the differences between the stable
    demographic vector (the dominant right eigenvector of \code{A})
    and the demographic distribution vector of the population (given
    by \code{vector}).\cr

    \code{Keyfitz.delta} will not work for reducible matrices and returns a warning
    for imprimitive matrices (although will not function for imprimitive matrices with
    nonzero imaginary components in the dominant eigenpair).

}
\value{
    Keyfitz's delta.
}
\references{
    Keyfitz (1968) Introduction to the Mathematics of Populations. Addison-Wesley.\cr
    Stott et al. (2010) Ecol. Lett., 14, 959-970.
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

    # Calculate Keyfitz's delta
    Keyfitz.delta(A, vector=initial)
}
\keyword{ Keyfitz }
\keyword{ distance }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }