\name{ inertia.tfsens }
\alias{ inertia.tfsens }
\alias{ inertia.tfsensmatrix }
\title{ Calculate sensitivity of inertia using transfer functions }
\description{
    Calculate the sensitivity of inertia to perturbations to a population
    projection matrix (PPM) model using a transfer function method.
}
\usage{
    inertia.tfsens(A, d=NULL, e=NULL, vector="n", bound=NULL,  
                   startval=0.001, tolerance=1e-10,
                   return.fit=FALSE, plot.fit=FALSE)
    inertia.tfsensmatrix(A, vector="n", bound=NULL,
                         startval=0.001, tolerance=1e-10)
}
\arguments{
  \item{A}{
    a square, primitive, nonnegative matrix of any dimension
  }
  \item{d,e}{
    numeric vectors that determine the perturbation structure (see details).
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate case-specific sensitvity of inertia.
  }
  \item{bound}{
    (optional) specifies whether an upper or lower bound should be calculated (see details).
  }
  \item{startval}{
    \code{inertia.tfsens} calculates the limit of the derivative of the transfer function as lambda approaches
    the dominant eigenvalue of \code{A} (see details).  \code{startval} provides a starting value for the
    algorithm: the smaller \code{startval} is, the quicker the algorithm should converge.
  }
  \item{tolerance}{
    the tolerance level for determining convergence (see details). 
  }
  \item{return.fit}{
    if \code{TRUE} (and only if \code{d} and \code{e} are specified), the lambda and sensitivity values obtained 
    in the convergence algorithm are returned alongside the calculated sensitivity of inertia.
  }
  \item{plot.fit}{
    if \code{TRUE} then convergence of the algorithm is plotted as sensitivity~lambda.
  }
}
\details{
    \code{inertia.tfsens} and \code{inertia.tfsensmatrix} differentiate a transfer function to find sensitivity 
    of inertia to perturbation.\cr

    If \code{vector="n"} then either \code{bound="upper"} or \code{bound="lower"}
    must be specified, so calculating the sensitivity of the upper or lower bound on population 
    inertia respectively.  Specifying \code{vector} overrides calculation of a bound, and will yield 
    sensitivity of case-specific inertia.\cr

    \code{inertia.tfsens} may be used to find sensitivity of a particular perturbation structure.  A desired 
    perturbation structure can be determined by \code{d\%*\%t(e)}.  Therefore, the rows to be perturbed
    are determined by \code{d} and the columns to be perturbed are determined by \code{e}.  The specific
    values in d and e will determine the relative perturbation magnitude.  So for example, if only entry
    [3,2] of a 3 by 3 matrix is to be perturbed, then \code{d = c(0,0,1)} and \code{e = c(0,1,0)}.  If entries 
    [3,2] and [3,3] are to be perturbed with the magnitude of perturbation to [3,2] half that of [3,3] then 
    \code{d = c(0,0,1)} and \code{e = c(0,0.5,1)}.  \code{d} and \code{e} may also be expressed as column
    vectors of class matrix, e.g. \code{d = matrix(c(0,0,1), ncol=1)}, \code{e = matrix(c(0,0.5,1), ncol=1)}.
    See Hodgson et al. (2006) for more information on perturbation structures.\cr

    \code{inertia.tfsensmatrix} returns a matrix of sensitivity values for observed transitions, where the sensitivity 
    of each matrix element is evaluated seperately.

    The formula used by \code{inertia.tfsens} and \code{inertia.tfsensmatrix} cannot be evaluated at lambda-max, therefore 
    it is necessary to find the limit of the formula as lambda approaches lambda-max.  This is done using a bisection method, 
    starting at a value of lambda-max + \code{startval}. \code{startval} should be small, to avoid the potential of false convergence.
    The algorithm continues until successive sensitivity calculations are within an accuracy of one another, determined by 
    \code{tolerance}: a \code{tolerance} of 1e-10 means that the sensitivity calculation should be accurate to 10 decimal places.
    However, as the limit approaches lambda-max, matrices become uninvertible (singular): if matrices are found to be singular
    then \code{tolerance} should be relaxed and made larger.

    For \code{inertia.tfsens}, there is an extra option to return and/or plot the above fitting process using \code{return.fit=TRUE} 
    and \code{plot.fit=TRUE} respectively.
}
\value{
    For \code{inertia.tfsens}, the sensitivity of inertia (or its bound) to the specified perturbation structure.  If \code{return.fit=TRUE},
    a list containing components:
  \item{sens}{
    the sensitivity of population inertia (or its bound) to the specified perturbation structure
  }
  \item{lambda.fit}{
    the lambda values obtained in the fitting process
  }
  \item{sens.fit}{
    the sensitivity values obtained in the fitting process
  }
    For \code{inertia.tfsensmatrix}, a matrix containing sensitivity of inertia to each seperate element of \code{A}. 
}
\references{
    Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Related:\cr
    \code{\link{inertia}}\cr

    Transfer function methods:\cr
    \code{\link{inertia.tfa}}, \code{\link{inertia.tfamatrix}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure    
    initial <- c(1,3,2)
    initial

    # Calculate the sensitivity matrix for the upper bound on inertia
    inertia.tfsensmatrix(A, bound="upper",tolerance=1e-7)

    # Calculate the sensitivity of simultaneous perturbation to 
    # A[1,3] and A[2,3] for the lower bound on inertia
    inertia.tfsens(A, d=c(1,0,0), e=c(0,1,1), bound="lower")

    # Calculate the sensitivity of simultaneous perturbation to 
    # A[1,3] and A[2,3] for specified initial stage structure
    # and return and plot the fitting process
    inertia.tfsens(A, d=c(1,0,0), e=c(0,1,1), vector=initial,tolerance=1e-7,
                   return.fit=TRUE,plot.fit=TRUE)
}
\keyword{ population }
\keyword{ inertia }
\keyword{ transfer function }
\keyword{ sensitivity }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }