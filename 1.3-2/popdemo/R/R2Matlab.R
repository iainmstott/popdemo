################################################################################
#' Convert matrices into Matlab style strings
#'
#' @description
#' Convert \R objects of class matrix into character strings that represent the 
#' matrix in a Matlab style
#'
#' @param A a numeric matrix of any dimension
#' @param noquote (optional) if \code{noquote=TRUE} then the returned character 
#' vector is printed without quotes.
#'
#' @details 
#' Matlab reads matrices using a unique one-line notation that can prove useful 
#' for storage in databases and importing multiple matrices into a program at 
#' once, amongst other applications.  This notation is by row, with "[" and "]" 
#' to specify the beginning and end of the matrix respectively, ";" to specify a 
#' new row and a space between each matrix element. Thus, the \R matrix created 
#' using \code{matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)} is 
#' equivalent to [0 1 2;0.5 0.1 0;0 0.6 0.6].
#' 
#' \code{R2Matlab} takes an \R object of class matrix converts it into a 
#' Matlab-style character string that may be useful for exporting into databases.
#'
#' @return 
#' Object of class character representing \code{A} in a Matlab style.
#'
#' @seealso
#' \code{\link{Matlab2R}}
#'
#' @examples
#'   # Create a 3x3 PPM
#'   ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )
#'
#'   # Code the matrix in a Matlab style
#'   R2Matlab(A)
#'
#'   # Print without quotes
#'   R2Matlab(A, noquote=TRUE)
#'
#' @concept Matlab
#'
#' @export R2Matlab
#'
R2Matlab <-
function(A,noquote=FALSE){
rows<-rep(0,nrow(A))
for(i in 1:nrow(A)){
    rows[i]<-paste(as.vector(A[i,]),collapse=" ")
}
matlabstr<-paste(noquote(rows),collapse=";")
matlabstr<-paste(noquote(c("[",matlabstr,"]")),collapse="")
if(noquote) matlabstr<-noquote(matlabstr)
return(matlabstr)
}

