#' @method print tfa
#' @export
print.tfa <- function(x,...){
    attr(x, "class")<-NULL
    print(x)
}

#' @method print tfam
#' @export
print.tfam <- function(x,...){
    attr(x,"class")<-NULL
    attr(x,"layout")<-NULL
    print(x)
}
