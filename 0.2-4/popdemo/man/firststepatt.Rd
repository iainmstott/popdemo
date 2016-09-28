\name{ firststepatt }
\alias{ firststepatt }
\title{ Calculate first-timestep attenuation }
\description{
    Calculate first-timestep attenuation for a specified 
    population projection matrix (PPM) model.
}
\usage{
    firststepatt(A, vector="n", return.N=FALSE)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate a case-specific first-timestep attenuation
  }
  \item{return.N}{
    (optional) if \code{TRUE}, returns population size in the first time 
    interval (including effects of asymptotic growth and initial population size), 
    alongside standardised first-timestep attenuation.
  }
}
\details{
    \code{firststepatt} returns a standardised measure of first-timestep 
    attenuation, so discounting the effects of both initial population size and asymoptotic 
    growth (Stott et al. 2011).\cr

    If \code{vector} is not specified then the bound on first-timestep attenuation
    (the smallest first-timestep attenuation that may be achieved) is returned, otherwise a 
    case-specific first-timestep attenuation for the specified matrix and demographic structure 
    is calculated.  Note that not all specified demographic structures will yield a 
    first-timestep attenuation: if the model does not attenuate in the first timestep then an 
    error is returned and \code{reactivity} should be used.\cr

    If \code{return.N=T} then the function also returns realised population size 
    (including the effects of asymptotic growth and initial population size).

    \code{firststepatt} works with imprimitive and irreducible matrices, but returns
    a warning in these cases.
}
\value{
    If \code{vector="n"}, the bound on first-timestep attenuation of \code{A}.\cr
    If \code{vector} is specified, the case-specific first-timestep attenuation of the model.\cr
    If \code{return.N=TRUE}, a list with components:\cr
  \item{firststepatt}{
    the bound on or case-specific first-timestep attenuation
  }
  \item{N}{
    the population size at the first timestep, including the effects 
    of initial population size and asymptotic growth.
  }
}
\references{
    Stott et al. (2011) Ecol. Lett., 14, 959-970.\cr
    Townley & Hodgson (2008) J. Appl. Ecol., 45, 1836-1839.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Other indices of transient density:\cr
    \code{\link{reactivity}}, \code{\link{maxamp}},
    \code{\link{maxatt}}, \code{\link{inertia}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(3,1,1)
    initial

    # Calculate the bound on first-timestep attenuation of A
    firststepatt(A)

    # Calculate case-specific first-timestep attenuation of A 
    # and initial
    firststepatt(A, vector=initial)

    # Calculate case-specific first-timestep attenuation of A 
    # and initial and return realised population size
    firststepatt(A, vector=initial, return.N=TRUE)
}
\keyword{ first-timestep attenuation }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }