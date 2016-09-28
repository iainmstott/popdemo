\name{ project }
\alias{ project }
\title{ Project population dynamics }
\description{
    Project dynamics of a specified population 
    projection matrix (PPM) model.
}
\usage{
    project(A, vector="n", time=100, standard.A=FALSE, 
            standard.vec=FALSE, return.vec=FALSE)
}
\arguments{
  \item{A}{
    a square, non-negative numeric matrix of any dimension.
  }
  \item{vector}{
    (optional) a specified initial age/stage distribution of class vector or 
    class matrix with which to project dynamics.
  }
  \item{time}{
    the number of projection intervals.
  }
  \item{standard.A}{
    (optional) if \code{TRUE}, standardises the PPM by scaling it by its dominant eigenvalue so that
    asymptotic dynamics of the projection are removed.
  }
  \item{standard.vec}{
    (optional) if \code{TRUE}, standardises \code{vector} to sum to 1 by scaling it by \code{sum(vector)}.
  }
  \item{return.vec}{
    (optional) if \code{TRUE}, returns the time series of vectors of demographic distribution
    as well as overall population size.
  }
}
\details{
    If \code{vector} is specified, then \code{project} will project the dynamics of the specified model using that vector.
    However, if \code{vector="n"}, then \code{project} will automatically project the set of stage-biased 
    vectors of \code{A}.  In practise, these projections are achieved using a set of standard basis vectors, 
    each with every element equal to 0, except for a single element that is equal to 1 (i.e. for a 3 by 3 
    matrix, the set of stage-biased vectors are: \code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}.)\cr

    Projections returned are of length \code{time+1}, as the first element represents the population
    at \code{t=0}.

    Projections have their own S3 plotting method to enable easy graphing.
}
\value{
    If \code{vector} is specified, a vector of population sizes of length \code{time+1}.\cr
    If \code{vector="n"}, a matrix of population projections: each column represents a single
    stage-biased projection and each column is of length \code{time+1}.\cr
    If \code{return.vec=TRUE}, a list with components:
  \item{N}{
    the vector or matrix of population sizes
  }
  \item{vec}{
    If \code{vector} is specified, a matrix of demographic vectors from projection of \code{vector} 
    through \code{A}.  Each column represents the densities of one life stage in the projection.\cr
    If \code{vector="n"}, an array of demographic vectors from projection of
    the set of stage-biased vectors through \code{A}.  The first dimension represents time (and is
    therefore equal to \code{time+1} in length).  The second dimension represents the densities of each
    stage (and is therefore equal to the dimension of \code{A} in length).  The third dimension
    represents each individual stage-biased projection (and is therefore also equal to the dimension of 
    \code{A} in length).  So for example, if we projected a 3 by 3 matrix for >10 time intervals (see examples)
    then the density of stage 3 in bias 2 at time 10 is found at element [11,3,2] (remembering that because 
    element 1 represents t=0, then t=10 is found at element 11); the time series of densities of stage 2
    in bias 1 is found using [,2,1]; the matrix of population vectors for bias 2 would be found using 
    [,,2].
  }    
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    S3 plotting method: \code{\link{plot.projection}}
}
\examples{
    # Create a 3x3 PPM
    A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
    A

    # Create an initial stage structure
    initial <- c(1,3,2)
    initial

    # Project stage-biased dynamics of A over 70 intervals
    pr <- project(A,time=70)
    pr
    plot(pr)

    # Select the projection of stage 2 bias
    pr[,2]

    # Project stage-biased dynamics of standardised A over 30 
    # intervals and return demographic vectors
    pr2 <- project(A, time=30, standard.A=TRUE, return.vec=TRUE)
    pr2
    plot(pr2)

    #Select the projection of stage 2 bias
    pr2$N[,2]

    # Select the density of stage 3 in bias 2 at time 10
    pr2$vec[11,3,2]

    # Select the time series of densities of stage 2 in bias 1
    pr2$vec[,2,1]

    #Select the matrix of population vectors for bias 2
    pr2$vec[,,2]

    # Project A over 50 intervals using a specified population structure
    pr3 <- project(A, vector=initial, time=50)
    pr3
    plot(pr3)

    # Project standardised dynamics of A over 10 intervals using 
    # standardised initial structure and return demographic vectors
    pr4 <- project(A, vector=initial, time=10, standard.vec=TRUE, 
                   standard.A=TRUE, return.vec=TRUE)
    pr4
    plot(pr4)
    
    # Select the time series for stage 1
    pr4$vec[,1]
    
}
\keyword{ project }
\keyword{ ecology }
\keyword{ demography }
\keyword{ population projection }
\keyword{ transient dynamics }
\keyword{ population projection matrix }
\keyword{ PPM }