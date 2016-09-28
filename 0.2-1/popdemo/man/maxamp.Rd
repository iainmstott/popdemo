\name{ maxamp }
\alias{ maxamp }
\title{ Calculate maximal amplification }
\description{
    Calculate maximal amplification for a specified 
    population projection matrix (PPM) model.
}
\usage{
    maxamp(A, vector = "n", return.N=FALSE, return.t=FALSE, return.stage=FALSE,
            conv.iterations=1e+5, conv.accuracy=1e-5)
}
\arguments{
  \item{A}{
    a square, primitive, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate a case-specific maximal amplification
  }
  \item{return.N}{
    (optional) if \code{TRUE}, returns population size at the point of maximal amplification 
    (including effects of asymptotic growth and initial population size), alongside standardised 
    maximal amplification.
  }
  \item{return.t}{
    (optional) if \code{TRUE}, returns the time at which maximal amplification occurs
    in the population projection.
  }
  \item{return.stage}{
    (optional) if \code{TRUE} and \code{vector="n"}, returns the stage that achieves the
    bound on maximal amplification
  }
  \item{conv.iterations}{
    the maximum number of iterations allowed when calulating convergence time (see details).
    Please see \code{iterations} in \code{\link{convergence.time}}.
  }
  \item{conv.accuracy}{
    the accuracy of convergence (see details).  Please see \code{accuracy} in 
    \code{\link{convergence.time}}.
  }
}
\details{
    \code{maxamp} returns a standardised measure of maximal amplification, so discounting 
    the effects of both initial population size and asymoptotic growth (Stott et al. 2011).  

    If \code{vector} is not specified then the bound on maximal amplification (the largest maximal 
    amplification that may be achieved) is returned, otherwise a case-specific maximal amplification 
    for the specified PPM and demographic structure is calculated.  Note that not all specified
    demographic structures will yield a maximal amplification: if the model does not amplify then an 
    error is returned.\cr

    Setting \code{return.N=T}, \code{return.t=T} and \code{return.stage=T} results in the function
    returning realised population size at maximal amplification (including the effects of asymptotic
    growth and initial population size), the time interval at which maximal amplification occurs and (if 
    \code{vector="n"}), the stage-bias that results in the bound on maximal amplification, respectively.
    NOTE that \code{N} is not indicative of maximum possible population size for a non-standardised model:
    merely the population size at the point of maximal amplification (i.e. largest 
    positive deviation from lambda-max).

    \code{max.amp} uses a simulation technique, using \code{\link{project}} to project the dynamics of 
    the model before evaluating maximum projected density over all t.  \code{conv.accuracy} and 
    \code{conv.iterations} are passed to \code{\link{convergence.time}}, which is used to find the point 
    of model convergence in order to ensure maximal amplification is correctly captured in model projection.\cr

    \code{maxamp} will not work for imprimitive or reducible matrices.
}
\value{
    If \code{vector="n"}, the bound on maximal amplification of \code{A}.\cr
    If \code{vector} is specified, the case-specific maximal amplification of the model.\cr
    If \code{return.N=TRUE}, \code{return.t=TRUE} and/or \code{return.stage=TRUE},
    a list with possible components:\cr
  \item{maxamp}{
    the bound on or case-specific maximal amplification
  }
  \item{N}{
    the population size at the point of maximal amplification, including the effects 
    of initial population size and asymptotic growth.  NOTE that \code{N} is not 
    indicative of maximum possible population size for a non-standardised model:
    merely the population size at the point of maximal amplification (i.e. largest 
    positive deviation from lambda-max).
  }
  \item{t}{
    the projection interval at which maximal amplification is achieved
  }
  \item{stage}{
    (only if \code{vector="n"}), the stage that achieves the bound on maximal amplification.
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
    \code{\link{reactivity}}, \code{\link{firststepatt}}
    \code{\link{maxatt}}, \code{\link{inertia}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Calculate the bound on maximal amplification of A
    maxamp(A)

    # Calculate the bound on maximal amplification of A and 
    # return the stage that achieves it
    maxamp(A, return.stage=TRUE)

    # Calculate case-specific maximal amplification of A
    # and initial
    maxamp(A, vector=initial)

    # Calculate case-specific maximal amplification of A
    # and initial and return realised population size and the 
    # time at which it is achieved
    maxamp(A, vector=initial, return.N=TRUE, return.t=TRUE)
}
\keyword{ maximal amplification }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }