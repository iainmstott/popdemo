\name{ inertia.tfa }
\alias{ inertia.tfa }
\title{ Analyse transfer function of population inertia }
\description{
    Analyse the transfer function of inertia of a given population 
    projection matrix (PPM) model and perturbation structure.    
}
\usage{
    inertia.tfa(A, d, e, vector="n", bound=NULL, prange, digits=1e-10)
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
    class matrix with which to calculate a case-specific inertia
  }
  \item{bound}{
    (optional) specifies whether an upper or lower bound should be calculated (see details).
  }
  \item{prange}{
    a numeric vector giving the continuous range of perturbation magnitude.
  }
  \item{digits}{
    specifies which values of lambda should be excluded to avoid a computationally singular system
    (see details).
  }
}
\details{
    \code{inertia.tfa} calculates the transfer function of inertia of a matrix given a perturbation 
    structure and a range of desired perturbation magnitude.  Currently \code{inertia.tfa} can 
    only work with rank-one, single-parameter perturbations (see Hodgson & Townley 2006).\cr

    If \code{vector="n"} then either \code{bound="upper"} or \code{bound="lower"}
    must be specified, so calculating the transfer function of the upper or lower bound on population 
    inertia respectively.  Specifying \code{vector} overrides calculation of a bound, and will yield a 
    transfer function of case-specific inertia.\cr

    The perturbation structure is determined by \code{d\%*\%t(e)}.  Therefore, the rows to be perturbed
    are determined by \code{d} and the columns to be perturbed are determined by \code{e}.  The specific
    values in d and e will determine the relative perturbation magnitude.  So for example, if only entry
    [3,2] of a 3 by 3 matrix is to be perturbed, then \code{d = c(0,0,1)} and \code{e = c(0,1,0)}.  If entries 
    [3,2] and [3,3] are to be perturbed with the magnitude of perturbation to [3,2] half that of [3,3] then 
    \code{d = c(0,0,1)} and \code{e = c(0,0.5,1)}.  \code{d} and \code{e} may also be expressed as column
    vectors of class matrix, e.g. \code{d = matrix(c(0,0,1), ncol=1)}, \code{e = matrix(c(0,0.5,1), ncol=1)}.
    See Hodgson et al. (2006) for more information on perturbation structures.\cr

    The perturbation magnitude is determined by \code{prange}, a numeric vector that gives the continuous
    range of perterbation magnitude to evaluate over.  This is usually a sequence (e.g. 
    \code{prange=seq(-0.1,0.1,0.001)}), but single transfer functions can be calculated using a single
    perturbation magnitude (e.g. \code{prange=-0.1}).  Because of the nature of the equation, used in calculating
    the transfer function, \code{prange} is used to find a range of lambda from which the perturbation magnitudes 
    are back-calculated and matched to their corresponding inertia, so that the output perturbation magnitude \code{p} 
    will match \code{prange} in length and range but not in numerical value (see value and Hodgson & Townley 2004 for
    more information).

    \code{inertia.tfa} uses the resolvent matrix in its calculation, which cannot be computed if any
    \code{lambdarange} are equal to the dominant eigenvalue of \code{A}.  \code{digits} specifies the
    values of lambda that should be excluded in order to avoid a computationally singular system.
    Any values equal to the dominant eigenvalue of \code{A} rounded to an accuracy of \code{digits}
    are excluded.  \code{digits} should only need to be changed when the system is found to be 
    computationally singular, in which case increasing \code{digits} should help to solve the problem.
    If not, then changing the range or step size of \code{prange} may help.\cr

    \code{inertia.tfa} will not work for imprimitive or reducible matrices.\cr

    There is an S3 plotting method available (see \code{\link{plot.tfa}} and examples below)
}
\value{
    A list containing numerical vectors:
  \item{p}{
    the perturbation magnitude.
  }
  \item{lambda}{
    the dominant eigenvalue of the perturbed matrix.
  }
  \item{inertia}{
    the inertia of the perturbed matrix model.
  }
    (Note that \code{p} will not usually be exactly the same as \code{prange} when specified, as the code
    calculates p for a given lambda rather than the other way around, with \code{prange} only used
    to determine max, min and number of lambda values to evaluate.)
}
\references{
    Townley & Hodgson (2004) J. Appl. Ecol., 41, 1155-1161.\cr
    Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    S3 plotting method:\cr
    \code{\link{plot.tfa}}\cr

    Related:\cr
    \code{\link{inertia}}, \code{\link{inertia.tfamatrix}}\cr

    Sensitivity methods:\cr
    \code{\link{inertia.tfsens}}, \code{\link{inertia.tfsensmatrix}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the transfer function of upper bound on inertia 
    # given a perturbation to A[3,2]
    transfer<-inertia.tfa(A, d=c(0,0,1), e=c(0,1,0), bound="upper",
                          prange=seq(-0.6,0.4,0.01))
    transfer

    # Plot the transfer function (defaults to p and inertia in
    # this case)
    plot(transfer)

    # Plot lambda and inertia
    plot(transfer, xvar="lambda", yvar="inertia")

    # Calculate the transfer function of case-specific inertia 
    # given perturbation to A[3,2] and A[3,3] with perturbation 
    # to A[3,2] half that of A[3,3]
    transfer2<-inertia.tfa(A, d=c(0,0,1), e=c(0,0.5,1), vector=initial,
                          prange=seq(-0.6,0.4,0.01))
    transfer2

    # Plot p and inertia
    plot(transfer2)

    # Plot lambda and inertia by hand
    plot(transfer$inertia~transfer$lambda,type="l", 
         xlab=expression(lambda),ylab="inertia")

}
\keyword{ population }
\keyword{ inertia }
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }