\name{ convergence.time }
\alias{ convergence.time }
\title{ Calculate time to convergence }
\description{
    Calculate the time to convergence for a
    specified population projection matrix (PPM) model.
}
\usage{
    convergence.time(A, vector = "n", accuracy=1e-2, iterations=1e+5)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate convergence time.
  }
  \item{accuracy}{
    The accuracy of convergence, expressed as a proportion (see details).
  }
  \item{iterations}{
    the maximum number of iterations of the model before the code breaks.  For slowly-converging 
    models and/or high specified convergence accuracy, this may need to be increased.
  }
}
\details{
    \code{convergence.time} works by simulating the model and manually determining when
    convergence to the given \code{accuracy} is reached.  Convergence is deemed to have been reached 
    when the growth of the model stays within the window determined by \code{accuracy} for 10*s 
    iterations of the model, with s equal to the dimension of \code{A}.  So for example, for an 8 by 8 
    matrix with dominant eigenvalue of 0.954 and convergence accuracy of 1e-2, growth of the model 
    must remain between 0.954*(1-1e-2)=0.9445 and 0.954*(1+1e-2)=0.9635 for 10*8=80 iterations of the 
    model for it to be considered 'converged'.\cr

    If \code{vector} is specified, then the convergence time of the projection of \code{vector}
    through \code{A} is returned.  However, if \code{vector="n"} then convergence times of the set of 
    stage-biased vectors is calculated.  In practise, these projections are achieved using a set of 
    standard basis vectors, each with every element equal to 0, except for a single element that is 
    equal to 1 (i.e. for a 3 by 3 matrix, the set of stage-biased vectors are: \code{c(1,0,0)}, 
    \code{c(0,1,0)} and \code{c(0,0,1)}.)\cr

    Due to the way in which convergence is defined, \code{convergence.time} can only properly
    work for strongly ergodic models.  Therefore, it will not function for imprimitive (therefore
    potentially weakly ergodic) or reducible (therefore potentially nonergodic) models.
}
\value{
    If \code{vector} is specified, the convergence time of the model.\cr
    If \code{vector} is not specified, a vector of convergence times for corresponding stage-biased
    projections (therefore, the length of the vector is equal to the dimension of \code{A}).
}
\references{
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

    # Calculate the convergence time of the 3 stage-biased 
    # populations within 0.1% of lambda-max
    convergence.time(A, accuracy=1e-3)

    # Calculate the convergence time of the projection of initial and A
    # to within 0.001% of lambda-max
    convergence.time(A, vector=initial, accuracy=1e-5)
}
\keyword{ convergence time }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }