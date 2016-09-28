\name{ tfa }
\alias{ tfa }
\title{ Analyse transfer function }
\description{
    Analyse the transfer function of a given population projection matrix 
    (PPM) and perturbation structure.    
}
\usage{
    tfa(A, d, e, prange=NULL, lambdarange=NULL, digits=1e-10)
}
\arguments{
  \item{A}{
    a square, irreducible, nonnegative matrix of any dimension
  }
  \item{d,e}{
    numeric vectors that determine the perturbation structure (see details).
  }
  \item{lambdarange}{
    a numeric vector giving the range of lambda values (asymptotic growth rates) to be achieved.
  }
  \item{prange}{
    a numeric vector giving the range of perturbation magnitude.
  }
  \item{digits}{
    specifies which values of lambda should be excluded to avoid a computationally singular system
    (see details).
  }
}
\details{
    \code{tfa} calculates the transfer function of a matrix given a perturbation structure and a range
    of desired lambda values or a range of desired perturbation magnitude.  Currently \code{tfa} can 
    only work with rank-one, single-parameter perturbations (see Hodgson & Townley 2004).\cr

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
    perturbation magnitude (e.g. \code{prange=-0.1}).  Because of the nature of the equation, \code{prange}
    is used to find a range of lambda from which the perturbation magnitudes are back-calculated and
    matched to their corresponding inertia, so that the output perturbation magnitude \code{p} will match
    \code{prange} in length and range but not in numerical value (see value).  Alternatively, a vector
    \code{lambdarange} can be specified, representing a range of desired lambda values from which the
    corresponding perturbation values will be calculated.  Only one of either \code{prange} or 
    \code{lambdarange} may be specified.  

    \code{tfa} uses the resolvent matrix in its calculation, which cannot be computed if any
    \code{lambdarange} are equal to the dominant eigenvalue of \code{A}.  \code{digits} specifies the
    values of lambda that should be excluded in order to avoid a computationally singular system.
    Any values equal to the dominant eigenvalue of \code{A} rounded to an accuracy of \code{digits}
    are excluded.  \code{digits} should only need to be changed when the system is found to be 
    computationally singular, in which case increasing \code{digits} should help to solve the problem.\cr

    \code{tfa} will not work for reducible matrices.\cr

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
    (Note that \code{p} will not be exactly the same as \code{prange} when specified, as the code
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
    \code{\link{tfamatrix}}\cr

    Sensitivity methods:\cr
    \code{\link{tfsens}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate the transfer function of A[3,2] given a range of lambda
    evals<-eigen(A)$values
    lmax<-which.max(Re(evals))
    lambda<-Re(evals[lmax])
    lambdarange <- seq(lambda-0.1, lambda+0.1, 0.01)
    transfer<-tfa(A, d=c(0,0,1), e=c(0,1,0), lambdarange=lambdarange)
    transfer

    # Plot the transfer function (defaults to lambda~p)
    plot(transfer)

    # Calculate the transfer function of perturbation to A[3,2] and A[3,3]
    # with perturbation to A[3,2] half that of A[3,3], given a range of 
    # perturbation values
    p<-seq(-0.6,0.4,0.01)
    transfer2<-tfa(A, d=c(0,0,1), e=c(0,0.5,1), prange=p)
    transfer2

    # Plot p and lambda by hand
    plot(transfer$lambda~transfer$p,type="l", xlab="p",ylab=expression(lambda))
}
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }