#' @export
print.projection <- function(x,...){
    attr(x,"bounds") <- NULL
    attr(x,"proj") <- NULL
    attr(x,"seq") <- NULL
    attr(x,"vec") <-NULL
    attr(x,"class")<-NULL
    print(x)
}

#' @export
print.tfa <- function(x,...){
    attr(x, "class")<-NULL
    print(x)
}

#' @export
print.tfam <- function(x,...){
    attr(x,"class")<-NULL
    attr(x,"layout")<-NULL
    print(x)
}
