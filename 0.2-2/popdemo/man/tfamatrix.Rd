\name{ tfamatrix }
\alias{ tfamatrix }
\title{ Analyse transfer function }
\description{
    Analyse the transfer function of a given population projection matrix 
    (PPM) and perturbation structure.    
}
\usage{
    tfamatrix(A, elementtype=NULL, Flim=c(-1,10), Plim=c(-1,10),
              plength=100, digits=1e-10)
}
\arguments{
  \item{A}{
    a square, irreducible, nonnegative matrix of any dimension
  }
  \item{elementtype}{
    (optional) a matrix describing the structure of \code{A}: \code{"P"} denotes elements bounded at 1,
    i.e. survival, growth, regression, \code{"F"} denotes elements not bounded at 1, i.e. fecundity,
    fission, and \code{NA} denotes absent elements (see details).
  }
  \item{Flim,Plim}{
    The perturbation ranges for \code{"F"} and \code{"P"} elements, expressed as a proportion of their
    magnitude  (see details).
  }
  \item{plength}{
    The desired length of the perturbation ranges.
  }
  \item{digits}{
    specifies which values of lambda should be excluded to avoid a computationally singular system
    (see details).
  }
}
\details{
    \code{tfamatrix} calculates an array of transfer functions for each nonzero element in \code{A}.
    Each element is perturbed individually.  The function is most useful for use with the S3 method
    \code{\link{plot.tfamatrix}} to visualise how perturbations affect the life cycle and 
    easily compare the effect on different transitions.\cr

    The structure of the perturbations are determined by \code{elementtype}, \code{Flim}, \code{Plim} and
    \code{plength}.  \code{elementtype} gives the nature of each transition, specifying whether perturbations
    should be bounded at 1 or not.  If \code{elementtype} is not directly specified, then the function
    assigns its own types, with those in the first row attributed \code{"F"}, and elsewhere in the 
    matrix being attributed \code{"F"} if the element >1 and \code{"P"} if the element is <=1.  \code{Flim} and
    \code{Plim} then determine the desired perturbation magnitude, expressed as a proportion of the magnitude of
    the elements, whilst plength determines the length of the perturbation vector.  So, if an "F" element is equal 
    to 0.5 and \code{Flim=c(-1,10)} and \code{plength=100} then the perturbation to that element is 
    \code{seq(-1*0.5,10*0.5,100-1) = seq(-0.5,5,100)}.  The process is the same for \code{"P"} elements, except that
    these are bounded to give a maximum perturbed value of 1 (as growth/survival elements cannot be greater
    than unity).  Both \code{"F"} and \code{"P"} elements are bounded to give a minimum perturbed value of 0.\cr

    \code{tfamatrix} uses \code{\link{tfa}} to calculate transfer functions.  \code{digits} is passed to 
    \code{tfa} to prevent the problem of singular matrices (see details in \code{\link{tfa}}).\cr

    \code{tfamatrix} will not work for reducible matrices.
}
\value{
    A list containing numerical arrays:
  \item{p}{
    the perturbation magnitude.
  }
  \item{lambda}{
    the dominant eigenvalue of the perturbed matrix.
  }
    The first and second dimensions of the arrays are equivalent to the first and second dimensions of the matrix.
    The third dimension of the arrays are therefore the vectors returned by \code{tfa}.  Therefore selecting e.g.
    $lambda[3,2,] will select the lambda values for element [3,2] of the matrix (for more information see examples).
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
    \code{\link{plot.tfamatrix}}\cr

    Related:\cr
    \code{\link{tfa}}\cr
    
    Sensitivity methods:\cr
    \code{\link{tfsens}}
}
\examples{
    # Create a 3x3 matrix
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Calculate the matrix of transfer functions using default arguments
    tfmat<-tfamatrix(A)
    tfmat

    # Plot the result
    plot(tfmat)
    
    # Plot the transfer function of element [3,2]
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
    plot(tfmat$lambda[3,2,]~tfmat$p[3,2,],xlab="p",ylab="lambda",type="l")

    # Create a new matrix with fission of adults
    B <- A
    B[2,3] <- 0.9
    B

    # Calculate the matrix of transfer functions using chosen arguments
    # that give the exact structure of the new matrix
    # and perturb a minimum of half the value of an element and
    # a maximum of double the value of an element
    etype <- matrix(c(NA, "F", "F", "P", "P", "F", NA, "P", "P"), 
                    ncol=3, byrow=TRUE)
    etype
    tfmat2 <- tfamatrix(B, elementtype=etype, Flim=c(-0.5,2),
                        Plim=c(-0.5,2))
    tfmat2

    # Plot the new matrix of transfer functions
    plot(tfmat2)
    
}
\keyword{ transfer function }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection matrix }
\keyword{ PPM }