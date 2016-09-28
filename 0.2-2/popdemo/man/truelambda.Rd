\name{ truelambda }
\alias{ truelambda }
\title{ Calculate simulated asymptotic growth }
\description{
    Calculate the true asymptotic growth of a specified population
    projection matrix (PPM) model using simulation.
}
\usage{
    truelambda(A, vector = "n", accuracy=1e-7, iterations=1e+5)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to project the model.
  }
  \item{accuracy}{
    the accuracy with which to determine convergence on asymptotic growth, expressed as a 
    proportion (see details).
  }
  \item{iterations}{
    the maximum number of iterations of the model before the code breaks.  For slowly-converging 
    models and/or high specified convergence accuracy, this may need to be increased.
  }
}
\details{
    \code{truelambda} works by simulating the given model and manually determining growth when
    convergence to the given \code{accuracy} is reached.  Convergence on an asymptotic growth 
    is deemed to have been reached when the growth of the model stays within the window determined by 
    \code{accuracy} for 10*s iterations of the model, with s equal to the dimension of \code{A}.  So 
    for example, projection of an 8 by 8 matrix with convergence accuracy of 1e-2 is deemed to have 
    converged on asymptotic growth when 10*8=80 consecutive iterations of the model have a growth 
    within 1-1e-2=0.99 (i.e. 99\%) and 1+1e-2=1.01 (i.e. 101\%) of each other.\cr
    
    If \code{vector} is specified, then the asymptotic growth of the projection of \code{vector}
    through \code{A} is returned.  If \code{vector="n"} then asymptotic growths of the set of stage-biased 
    vectors is calculated.  In practise, these projections are achieved using a set of standard basis vectors,
    each with every element equal to 0, except for a single element that is equal to 1 (i.e. for a 3 by 3 
    matrix, the set of stage-biased vectors are: \code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}.)
}
\value{
    If \code{vector} is specified, a numeric vector of length 2 giving the range in which asymptotic
    growth of the model lies.\cr
    
    If \code{vector} is not specified, a 2-column matrix with each row giving the range in which
    asymptotic growth lies for its corresponding stage-biased projection.  Therefore, the number of 
    rows is equal to the dimension of \code{A}.
}
\references{
    Stott et al. (2010) Methods Ecol. Evol., 1, 242-252.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\examples{
    # Create a 3x3 irreducible PPM
    A <- matrix(c(0,0,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the true asymptotic growth of the stage-biased
    # projections of A
    truelambda(A)

    # Calculate the true asymptotic growth of the projection of
    # A and initial
    truelambda(A, vector=initial)

    # Create a 3x3 reducible, nonergodic PPM
    B<-A
    B[3,2] <- 0
    B

    # Calculate the true asymptotic growth of the 3 stage-biased
    # projections of B
    truelambda(B)
}
\keyword{ asymptotic growth }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }