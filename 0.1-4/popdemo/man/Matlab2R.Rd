\name{ Matlab2R }
\alias{ Matlab2R }
\title{ Read Matlab style matrices into R }
\description{
    Read a matrix coded in a Matlab style into \R to
    create an object of class matrix
}
\usage{
    Matlab2R(M)
}
\arguments{
  \item{M}{
    an object of class character that represents a numeric matrix coded in a Matlab style.
  }
}
\details{
    Matlab reads matrices using a unique one-line notation that can prove useful for storage in databases
    and importing multiple matrices into a program at once, amongst other applications.  This notation 
    is by row, with "[" and "]" to specify the beginning and end of the matrix respectively, ";" to 
    specify a new row and a space between each matrix element.  Thus, the \R matrix created using 
    \code{matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)} may is equivalent to
    [0 1 2;0.5 0.1 0;0 0.6 0.6].\cr

    \code{Matlab2R} takes a Matlab-coded matrix expressed as a character string and uses string-splitting 
    and string-building techniques to convert it into an \R object of class matrix.  As well as providing a 
    simpler means of matrix notation in \R, it also enables simultaneous import of multiple matrices of 
    varying dimensions, using comma-seperated dataframes and tables.
}
\value{
    An object of class matrix.
}
\author{
    Stott, I., Hodgson, D. J., Townley, S.
}
\seealso{
    \code{\link{R2Matlab}}
}
\examples{
    # Create a 3x3 PPM using Matlab2R
    A<-Matlab2R("[0 1 2;0.5 0.1 0;0 0.6 0.6]")
    A
}
\keyword{ Matlab2R }
\keyword{ Matlab }
\keyword{ matrix }