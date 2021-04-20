################################################################################
#' Calculate asymptotic growth
#'
#' @description
#' Dominant eigenstuff of a population matrix projection model.
#'
#' @param A a square, nonnegative numeric matrix of any dimension.
#' @param what what components of the dominant eigenstuff should be returned. 
#' A character vector, which may include:
#' \describe{
#'  \item{\code{"lambda"}}{the dominant eigenvalue (lambda)}
#'  \item{\code{"ss"}}{the dominant right eigenvector (stable stage)}
#'  \item{\code{"rv"}}{the dominant left eigenvector (reproductive value)}
#' }
#' the default, \code{"all"}, returns all of the above.
#' @param check (logical) determines whether the dominant eigenvalue is 
#' checked for nonzero imaginary component, and largest absolute value. If
#' either of these occur, then the dominant eigenvalue may not be described as 
#' truly dominant.
#'
#' @details 
#' \code{eigs} gives the dominant eigenstuff of a population projection model. 
#' This includes the dominant eigenvalue (asymptotic population growth), the 
#' dominant right eigenvector (stable age/stage distribution), and the dominant
#' left eigenvector (reproductive value). The dominant eigenvalue is the 
#' eigenvalue with the largest real component, and the dominant eigenvectors are 
#' the eigenvectors that correspond to this eigenvalue. If the matrix is 
#' reducible, then there may be other real or complex eigenvalues whose absolute
#' value are equal in magnitude to that of the dominant eigenvalue. In this case,
#' \code{eigs} returns the first one, and gives a warning "More than one eigenvalues 
#' have equal absolute magnitude", for information.
#' 
#' @return 
#' A list with possible components that depends on the contents of \code{what}:
#' \describe{
#' \item{lambda}{ 
#' the dominant eigenvalue, which describes asymptotic population growth (if A 
#' is primitive; see \code{\link{isPrimitive}}). A real, nonnegative numeric 
#' vector of length 1. 
#' }
#' \item{ss}{
#' the dominant right eigenvector, which describes the stable age/stage structure
#' (if \code{A} is primitive; see \code{\link{isPrimitive}}). A real, nonnegative 
#' numeric vector equal to the dimension of \code{A} in length, scaled to sum to 1.
#' }
#' \item{rv}{
#' the dominant left eigenvector, which describes the reproductive value (if 
#' \code{A} is primitive; see \code{\link{isPrimitive}}). A real, nonnegative 
#' numeric vector equal to the dimension of \code{A} in length, scaled so that 
#' rv%*%ss equals 1.
#' }
#' }
#' 
#' If only one of these components is returned, then the value is not a list, but 
#' a single numeric vector.
#'
#' @family Eigenstuff 
#' @family GrowthMeasures
#'
#' @examples
#'   # load the desert tortoise data
#'   data(Tort)
#'   
#'   # find the dominant eigenvalue
#'   eigs(Tort, "lambda")
#'   
#'   #find the stable stage structure
#'   eigs(Tort, "ss")
#'   
#'   #find the reproductive value
#'   eigs(Tort, "rv")
#'   
#'   #find both dominant eigenvectors
#'   eigs(Tort, c("ss","rv"))
#'   
#'   #find all eigenstuff
#'   eigs(Tort)
#'
#' @concept 
#' asymptotic growth population structure reproductive value eigenvalues
#'
#' @export eigs
#'
eigs <-
function(A, what = "all", check = TRUE){
ifelse("lambda" %in% what, val <- TRUE, val <- FALSE)
ifelse("ss" %in% what, rvec <- TRUE, rvec <- FALSE)
ifelse("rv" %in% what, lvec <- TRUE, lvec <- FALSE)
if("all" %in% what){
    val <- TRUE
    rvec <- TRUE
    lvec <- TRUE
}
if(!val & !rvec & !lvec) stop('"what" does not contain the right information')
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
eigstuff <- eigen(A)
lmax <- which.max(Re(eigstuff$values))
lambda <- Re(eigstuff$values[lmax])
if(check){
    if(any(abs(eigstuff$values[-lmax]) >= lambda)){
        warning("More than one eigenvalues have equal absolute magnitude")
    }
    if(Im(eigstuff$values[lmax]) > 0){
        warning("'Dominant' eigenvalue has nonzero imaginary component")
    }
}
if(any(rvec, lvec)){
    ss <- Re(eigstuff$vectors[,lmax])
    ss <- ss / sum(ss)
}
if(lvec){
    teigstuff <- eigen(t(A))
    rv <- Re(teigstuff$vectors[,lmax])
    rvss <- as.numeric(rv%*%ss)
    rv <- rv / rvss
}
if(val&!rvec&!lvec) return(lambda)
if(!val&rvec&!lvec) return(ss)
if(!val&!rvec&lvec) return(rv)
if(val&rvec&!lvec) return(list(lambda = lambda, ss = ss))
if(val&!rvec&lvec) return(list(lambda = lambda, rv = rv))
if(!val&rvec&lvec) return(list(ss = ss, rv = rv))
if(val&rvec&lvec) return(list(lambda = lambda, ss = ss, rv = rv))
}

