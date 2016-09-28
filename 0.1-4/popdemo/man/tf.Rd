\name{ tf }
\alias{ tf }
\title{ Calculate transfer function }
\description{
    Calculate the transfer function of a given matrix and perturbation structure.    
}
\usage{
    tf(A, z, b, c, exp=-1)
}
\arguments{
  \item{A}{
    a matrix of any dimension.
  }
  \item{z}{
    a numeric vector representing the set of real values over which to analyse the transfer function.
  }
  \item{b,c}{
    a numeric column vectors of class matrix that determines the perturbation structure (see details).
  }
  \item{exp}{
    the exponent to raise the resolvent to.  Normally, \code{exp=-1} but for some usages an exponent 
    of larger magnitude may be needed.
  }
}
\details{
    \code{tf} is primarily designed to be called by \code{\link{tfa}}, \code{\link{tfsens}},
    \code{\link{inertia.tfa}} and \code{\link{inertia.tfsens}} for evaluating transfer functions and 
    sensitivities of population projection matrix models.  However, it can also function as a 
    standalone means of calculating a transfer function.  Currently \code{tf} can only work with 
    rank-one, single-parameter perturbations.\cr

    For \code{tf} to work properly, the dimensions of \code{A}, \code{b} and \code{c} must match.
    If \code{A} is an n by m matrix, then \code{b} must be a column vector of length m and \code{c} 
    must be a column vector of length n.\cr

    The perturbation structure is determined by \code{b\%*\%t(c)}.  Therefore, the rows to be perturbed
    are determined by \code{b} and the columns to be perturbed are determined by \code{c}.  The specific
    values in \code{b} and \code{c} will determine the relative perturbation magnitude.  So for example, 
    if only entry [3,2] of a 3 by 3 matrix is to be perturbed, then \code{b = matrix(c(0,0,1), ncol=1)}
    and \code{c = matrix(c(0,1,0), ncol=1)}.  If entries [3,2] and [3,3] are to be perturbed with the 
    magnitude of perturbation to [3,2] half that of [3,3] then \code{b = matrix(c(0,0,1), ncol=1)} and 
    \code{c = matrix(c(0,0.5,1), ncol=1)}.  See Hodgson et al. (2006) for more information on perturbation 
    structures.\cr
}
\value{
    a vector of perturbation magnitudes whose INVERSE will achieve the dominant eigenvalues given by \code{z}
    for the triple \code{A}, \code{b}, \code{c}.
}
\references{
    Townley & Hodgson (2004) J. Appl. Ecol., 41, 1155-1161.\cr
    Hodgson et al. (2006) J. Theor. Biol., 70, 214-224.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{tfa}}, \code{\link{tfsens}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate the transfer function given a range of z and a perturbation structure
    evals<-eigen(A)$values
    lmax<-which.max(Re(evals))
    lambda<-Re(evals[lmax])
    z <- seq(lambda-0.04, lambda+0.04, 0.01)
    z <- z[ z != lambda]
    b <- matrix(c(0,0,1), ncol=1)
    c <- matrix(c(0,0,1), ncol=1)
    p<-1/tf(A, z=z, b=b, c=c)
    p

    ## Plot p and z
    plot(z~p,type="l", xlab="p",ylab=expression(lambda[max]))
}
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }