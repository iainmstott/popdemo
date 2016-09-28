\name{ tfsens }
\alias{ tfsens }
\alias{ tfsensmatrix }
\title{ Calculate sensitivity using transfer functions }
\description{
    Calculate the sensitivity of perturbations to a population
    projection matrix using a transfer function method.
}
\usage{
    tfsens(A, d=NULL, e=NULL, startval=0.001, tolerance=1e-10, 
           return.fit=FALSE, plot.fit=FALSE)
    tfsensmatrix(A, startval=0.001, tolerance=1e-10)
}
\arguments{
  \item{A}{
    a square, nonnegative matrix of any dimension.
  }
  \item{d,e}{
    numeric vectors that determine the perturbation structure (see details).
  }
  \item{startval}{
    \code{tfsens} calculates the limit of the derivative of the transfer function as lambda approaches
    the dominant eigenvalue of \code{A} (see details).  \code{startval} provides a starting value for the
    algorithm: the smaller \code{startval} is, the quicker the algorithm should converge.
  }
  \item{tolerance}{
    the tolerance level for determining convergence (see details). 
  }
  \item{return.fit}{
    if \code{TRUE} the lambda and sensitivity values obtained in the convergence algorithm are returned 
    alongside the calculated sensitivity of inertia.
  }
  \item{plot.fit}{
    if \code{TRUE} then convergence of the algorithm is plotted as sensitivity~lambda.
  }
}
\details{
    \code{tfsens} and \code{tfsensmatrix} differentiate a transfer function to find sensitivity to perturbation.  
    This provides an alternative method to using matrix eigenvectors to calculate the sensitivity
    matrix and is useful as it may utilise more complex perturbation structures.\cr

    \code{tfsens} may be used to find sensitivity of a particular perturbation structure. A desired 
    perturbation structure is determined by \code{d\%*\%t(e)}.  Therefore, the rows to be perturbed
    are determined by \code{d} and the columns to be perturbed are determined by \code{e}.  The specific
    values in d and e will determine the relative perturbation magnitude.  So for example, if only entry
    [3,2] of a 3 by 3 matrix is to be perturbed, then \code{d = c(0,0,1)} and \code{e = c(0,1,0)}.  If entries 
    [3,2] and [3,3] are to be perturbed with the magnitude of perturbation to [3,2] half that of [3,3] then 
    \code{d = c(0,0,1)} and \code{e = c(0,0.5,1)}.  \code{d} and \code{e} may also be expressed as column
    vectors of class matrix, e.g. \code{d = matrix(c(0,0,1), ncol=1)}, \code{e = matrix(c(0,0.5,1), ncol=1)}.
    See Hodgson et al. (2006) for more information on perturbation structures.\cr

    \code{tfsensmatrix} returns a matrix of sensitivity values for observed transitions (similar to that obtained when using
    \code{\link{sens}} to evaluate sensitivity using eigenvectors), where the sensitivity of each matrix element
    is evaluated seperately.

    The formula used by \code{tfsens} and \code{tfsensmatrix} cannot be evaluated at lambda-max, therefore it is necessary 
    to find the limit of the formula as lambda approaches lambda-max.  This is done using a bisection method, starting at a 
    value of lambda-max + \code{startval}. \code{startval} should be small, to avoid the potential of false convergence.
    The algorithm continues until successive sensitivity calculations are within an accuracy of one another, determined by 
    \code{tolerance}: a \code{tolerance} of 1e-10 means that the sensitivity calculation should be accurate to 10 decimal places.
    However, as the limit approaches lambda-max, matrices become uninvertible (singular): if matrices are found to be singular
    then \code{tolerance} should be relaxed and made larger.  

    For \code{tfsens}, there is an extra option to return and/or plot the above fitting process using \code{return.fit=TRUE} 
    and \code{plot.fit=TRUE} respectively.
}
\value{
    For \code{tfsens}, the sensitivity of lambda-max to the specified perturbation structure.  If \code{return.fit=TRUE}
    a list containing components:\cr
  \item{sens}{
    the sensitivity of lambda-max to the specified perturbation structure
  }
  \item{lambda.fit}{
    the lambda values obtained in the fitting process
  }
  \item{sens.fit}{
    the sensitivity values obtained in the fitting process
  }
    For \code{tfsensmatrix}, a matrix containing sensitivity of lambda-max to each seperate element of \code{A}.
}
\references{
    Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Transfer function methods:\cr
    \code{\link{tfa}}, \code{\link{tfamatrix}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate the sensitivity matrix
    tfsensmatrix(A)

    # Calculate the sensitivity of simultaneous perturbation to 
    # A[1,2] and A[1,3]
    tfsens(A, d=c(1,0,0), e=c(0,1,1))

    # Calculate the sensitivity of simultaneous perturbation to 
    # A[1,2] and A[1,3] and return and plot the fitting process
    tfsens(A, d=c(1,0,0), e=c(0,1,1),
           return.fit=TRUE, plot.fit=TRUE)
}
\keyword{ transfer function }
\keyword{ sensitivity }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }