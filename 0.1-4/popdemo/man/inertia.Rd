\name{ inertia }
\alias{ inertia }
\title{ Calculate population inertia }
\description{
    Calculate population inertia for a specified 
    population projection matrix (PPM) model.
}
\usage{
    inertia(A, vector = "n", bound=NULL, return.N = FALSE, t=NULL)
}
\arguments{
  \item{A}{
    a square, primitive, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate a case-specific inertia
  }
  \item{bound}{
    (optional) specifies whether an upper or lower bound should be calculated (see details).
  }
  \item{return.N}{
    (optional) if \code{TRUE}, returns population size for specified \code{t} (including effects 
    of asymptotic growth and initial population size), alongside standardised inertia.
  }
  \item{t}{
    (optional) the projection interval at which \code{N} is to be calculated.  Calculation of \code{N} 
    is only accurate for \code{t} where the model has converged.
  }
}
\details{
    \code{inertia} returns a standardised measure of population inertia (Koons et al. 2007), so 
    discounting the effects of both initial population size and asymoptotic growth (Stott et al.
    2011).\cr

    If \code{vector="n"} then either \code{bound="upper"} or \code{bound="lower"}
    must be specified, so calculating the upper or lower bound on population inertia respectively
    (i.e. the largest and smallest values that inertia may take).  Specifying \code{vector} overrides 
    calculation of a bound, and will yield a case-specific value for inertia.\cr

    \code{inertia} will not work with imprimitive or reducible matrices.
}
\value{
    If \code{vector="n"}, the upper bound on inertia of \code{A} if \code{bound="upper"}
    and the lower bound on inertia of \code{A} if \code{bound="lower"}.\cr
    If \code{vector} is specified, the case-specific inertia of the model.\cr
    If \code{return.N=TRUE} and \code{t} is specified, a list with components:
  \item{inertia}{
    the bound on or case-specific inertia
  }
  \item{N}{
    the population size at specified \code{t}.
  }
}
\references{
    Koons et al. (2007) Ecology, 88, 2867-2867.\cr
    Stott et al. (2011) Ecol. Lett., 14, 959-970.\cr
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Other indices of transient density:\cr
    \code{\link{reactivity}}, \code{\link{firststepatt}},
    \code{\link{maxamp}}, \code{\link{maxatt}}\cr

    Transfer function methods:\cr
    \code{\link{inertia.tfa}}, \code{\link{inertia.tfamatrix}}\cr

    Sensitivity methods:
    \code{\link{inertia.tfsens}}, \code{\link{inertia.tfsensmatrix}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the upper bound on inertia of A
    inertia(A,bound="upper")

    # Calculate the lower bound on inertia of A
    inertia(A,bound="lower")

    # Calculate case-specific inertia of A and initial
    inertia(A, vector=initial)

    # Calculate case-specific inertia of A and initial and 
    # return realised population size at t=25
    inertia(A, vector=initial, return.N=TRUE, t=25)
}
\keyword{ inertia }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }