\name{ dr }
\alias{ dr }
\title{ Calculate damping ratio }
\description{
    Calculate the damping ratio of a given population
    projection matrix (PPM).
}
\usage{
    dr(A, return.time=FALSE, x=10)
}
\arguments{
  \item{A}{
    a square, irreducible, non-negative numeric matrix of any dimension
  }
  \item{return.time}{
    (optional) a logical argument determining whether an estimated
    convergence time should be returned.
  }
  \item{x}{
    (optional) the logarithm used in determining estimated time to convergence (see details).
  }
}
\details{
    The damping ratio is calculated as the ratio of the dominant eigenvalue to the
    modulus of the largest subdominant eigenvalue.  Time to convergence can be found
    by calculating \code{log(dr)/log(x)}, which is the time taken for the dominant eigenvalue
    to become \code{x} times larger than the largest subdominant eigenvalue.
}
\value{
    If \code{return.time=FALSE}, the damping ratio of \code{A}.\cr
    If \code{return.time=TRUE}, a list containing components:
  \item{dr}{
    the damping ratio of \code{A}
  }
  \item{t}{
    the estimated time to convergence.
  }
}
\references{
    Caswell (2001) Matrix Population Models 2nd. ed. Sinauer.\cr
    Stott et al. (2010) Ecol. Lett., 14, 959-970.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate damping ratio
    dr(A)

    # Calculate damping ratio and time to convergence using a 
    # multiple of 10
    dr(A, return.time=TRUE, x=10)
}
\keyword{ damping ratio }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }