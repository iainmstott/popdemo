\name{ reactivity }
\alias{ reactivity }
\title{ Calculate reactivity }
\description{
    Calculate reactivity (first-timestep amplification) for a specified 
    population projection matrix (PPM) model.
}
\usage{
    reactivity(A, vector = "n", return.N = FALSE)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate a case-specific reactivity
  }
  \item{return.N}{
    (optional) if \code{TRUE}, returns population size in the first time 
    interval (including effects of asymptotic growth and initial population size), 
    alongside standardised reactivity.
  }
}
\details{
    \code{reactivity} is a standardised measure of first-timestep amplification, so 
    discounts the effects of both initial population size and asymoptotic growth (Stott
    et al. 2011)\cr

    If \code{vector} is not specified then the bound on reactivity (the largest 
    reactivity that may be achieved) is returned, otherwise a case-specific reactivity 
    for the specified PPM and demographic structure is calculated.  Note that
    not all specified demographic structures will yield a reactivity: if the 
    model does not amplify in the first timestep then an error is returned and 
    \code{\link{firststepatt}} should be used.\cr

    If \code{return.N=T} then the function also returns realised population size 
    (including the effects of asymptotic growth and initial population size).

    \code{reactivity} works with imprimitive and irreducible matrices, but returns
    a warning in these cases.
}
\value{
    If \code{vector="n"}, the bound on reactivity of \code{A}.\cr
    If \code{vector} is specified, the case-specific reactivity of the model.\cr
    If \code{return.N=TRUE}, a list with components:\cr
  \item{reactivity}{
    the bound on or case-specific reactivity
  }
  \item{N}{
    the population size at the first timestep, including the effects 
    of initial population size and asymptotic growth.
  }
}
\references{
    Neubert & Caswell (1997) Ecology, 78, 653-665.\cr   
    Stott et al. (2011) Ecol. Lett., 14, 959-970.\cr
    Townley & Hodgson (2008) J. Appl. Ecol., 45, 1836-1839.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Other indices of transient density:\cr
    \code{\link{firststepatt}}, \code{\link{maxamp}}
    \code{\link{maxatt}}, \code{\link{inertia}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the bound on reactivity of A
    reactivity(A)

    # Calculate case-specific reactivity of A
    # when projected using specific demographic structure
    reactivity(A, vector=initial)

    # Calculate case-specific reactivity of A
    # and initial and return realised population size
    reactivity(A, vector=initial, return.N=TRUE)
}
\keyword{ reactivity }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }