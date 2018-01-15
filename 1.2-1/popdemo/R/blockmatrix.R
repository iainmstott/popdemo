################################################################################
#' Block-permute a reducible matrix
#'
#' @description
#' Conjugate a reducible matrix into block upper triangular form
#'
#' @param A a square, reducible, non-negative numeric matrix of any dimension
#'
#' @details 
#' Any square, reducible, non-negative matrix may have its rows and columns 
#' conjugated so that it takes a block upper triangular structure, with 
#' irreducible square submatrices on the diagonal, zero submatrices in the 
#' lower triangle and non-negative submatrices in the upper triangle (Caswell 
#' 2001; Stott et al. 2010). \code{blockmatrix} permutes the rows and columns 
#' of a reducible matrix into this form, which enables further evaluation (e.g. 
#' computation of eigenvalues of submatrices).
#'
#' @return 
#' a list containing components:
#' \describe{
#' \item{\code{blockmatrix}}{ the block-permuted matrix. }
#' \item{\code{stage.order}}{ the permutation of rows/columns of \code{A} in the
#' block-permuted matrix. }
#' }
#'
#' @references
#' Caswell (2001) Matrix population models 2nd ed. Sinauer.\cr
#' Stott et al. (2010) Methods. Ecol. Evol., 1, 242-252.
#'
#' @examples
#'   # Create a 3x3 reducible PPM
#'   A <- matrix(c(0,1,0,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3)
#'   dimnames(A) <- list(c("Juv", "Pre-R", "R"), c("Juv", "Pre-R", "R"))
#'   A
#'
#'   # Block-permute the matrix
#'   blockmatrix(A)
#' 
#' @concept
#' reducible irreducible submatrix permutation conjugation permute conjugate 
#' arrange rearrange row column
#'
#' @export blockmatrix
#' @importFrom expm "%^%"
#'
#
#Can be vectorised
#
blockmatrix <-
function(A){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(isIrreducible(A)) stop("Cannot compute block-permutation for irreducible matrix")
stage.order<-1:order
powermatrix<-(diag(1,order)+A)%^%(order-1)
zeroes1<-numeric(order)
for(i in 1:order){
    zeroes1[i]<-length(which(powermatrix[i,]==0))
}
stage.order1<-order(zeroes1)
blockpowermatrix<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        blockpowermatrix[i,j]<-powermatrix[stage.order1[i],stage.order1[j]]
    }
}
zeroes2<-numeric(order)
for(i in 1:order){
    if(blockpowermatrix[i,1]==0){
        while(blockpowermatrix[i,1+zeroes2[i]]==0){
            zeroes2[i]<-zeroes2[i]+1
        }
    }
}
stage.order2<-order(zeroes2)
stage.order<-stage.order[stage.order1][stage.order2]
blockmatrix<-matrix(0,order,order)
for(i in 1:order){
    for(j in 1:order){
        blockmatrix[i,j]<-A[stage.order[i],stage.order[j]]
    }
}
dimnames(blockmatrix)<-lapply(dimnames(A),function(x){x[stage.order]})
return(list(blockmatrix=blockmatrix,stage.order=stage.order))
}
