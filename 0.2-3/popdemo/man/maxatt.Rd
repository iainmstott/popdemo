\name{ maxatt }
\alias{ maxatt }
\title{ Calculate maximal attenuation }
\description{
    Calculate maximal attenuation for a specified 
    population projection matrix (PPM) model.
}
\usage{
    maxatt(A, vector = "n", return.N = FALSE,return.t=FALSE, return.stage=FALSE,
            conv.iterations=1e+5, conv.accuracy=1e-5)
}
\arguments{
  \item{A}{
    a square, primitive, non-negative numeric matrix of any dimension
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to calculate a case-specific maximal attenuation
  }
  \item{return.N}{
    (optional) if \code{TRUE}, returns population size at the point of maximal attenuation 
    (including effects of asymptotic growth and initial population size), alongside standardised 
    maximal attenuation.
  }
  \item{return.t}{
    (optional) if \code{TRUE}, returns the time at which maximal attenuation occurs
    in the population projection.
  }
  \item{return.stage}{
    (optional) if \code{TRUE} and \code{vector="n"}, returns the stage that achieves the
    bound on maximal attenuation
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
    \code{maxatt} returns a standardised measure of maximal attenuation, so discounting 
    the effects of both initial population size and asymoptotic growth (Stott et al. 2011).  

    If \code{vector} is not specified then the bound on maximal attenuation (the largest maximal 
    attenuation that may be achieved) is returned, otherwise a case-specific maximal attenuation 
    for the specified PPM and demographic structure is calculated.  Note that not all specified
    demographic structures will yield a maximal attenuation: if the model does not attenuate then an 
    error is returned.\cr

    Setting \code{return.N=T}, \code{return.t=T} and \code{return.stage=T} results in the function
    returning realised population size at maximal attenuation (including the effects of asymptotic
    growth and initial population size), the time interval at which maximal attenuation occurs and (if 
    \code{vector="n"}), the stage-bias that results in the bound on maximal attenuation, respectively.
    NOTE that \code{N} is not indicative of minimum possible population size for a non-standardised model:
    merely the population size at the point of maximal attenuation (i.e. largest 
    negative deviation from lambda-max).

    \code{max.att} uses a simulation technique, using \code{\link{project}} to project the dynamics of 
    the model before evaluating minimum projected density over all t.  \code{conv.accuracy} and 
    \code{conv.iterations} are passed to \code{\link{convergence.time}}, which is used to find the point 
    of model convergence in order to ensure maximal attenuation is correctly captured in model projection.\cr

    \code{maxatt} will not work for imprimitive or reducible matrices.
}
\value{
    If \code{vector="n"}, the bound on maximal attenuation of \code{A}.\cr
    If \code{vector} is specified, the case-specific maximal attenuation of the model.\cr
    In addition, if \code{return.N=TRUE}, \code{return.t=TRUE} and/or \code{return.stage=TRUE},
    a list with possible components:\cr
  \item{maxamp}{
    the bound on or case-specific maximal attenuation
  }
  \item{N}{
    the population size at the point of maximal attenuation, including the effects 
    of initial population size and asymptotic growth.  NOTE that \code{N} is not 
    indicative of minimum possible population size for a non-standardised model:
    merely the population size at the point of maximal attenuation (i.e. largest 
    negative deviation from lambda-max).
  }
  \item{t}{
    the projection interval at which maximal attenuation is achieved
  }
  \item{stage}{
    (only if \code{vector="n"}), the stage that achieves the bound on maximal attenuation.
  }
}
\references{
    Stott et al. (2011) Ecol. Lett., 14, 959-970.\cr
    Townley & Hodgson (2008) J. Appl. Ecol., 45, 1836-1839.  
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    Other indices of transient density:\cr
    \code{\link{reactivity}}, \code{\link{firststepatt}},
    \code{\link{maxamp}}, \code{\link{inertia}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(3,1,1)
    initial

    # Calculate the bound on maximal attenuation of A
    maxatt(A)

    # Calculate the bound on maximal attenuation of A and 
    # return the stage that achieves it
    maxatt(A,return.stage=TRUE)

    # Calculate case-specific maximal attenuation of A
    # and initial
    maxatt(A, vector=initial)

    # Calculate case-specific maximal attenuation of A
    # and initial and return realised population size and the 
    # time at which it is achieved
    maxatt(A, vector=initial, return.N=TRUE, return.t=TRUE)
}
\keyword{ maximal attenuation }
\keyword{ ecology }
\keyword{ demography }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }