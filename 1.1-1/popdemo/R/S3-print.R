#' @export
print.projection <- function(x,...){
    attr(x,"class")<-NULL
    attr(x,"bounds") <- NULL
    attr(x,"vec") <-NULL
    print(x)
}

#' @export
print.tfa <- function(x,...){
    attributes(x)$class<-NULL
    print(x)
}

#' @export
print.tfam <- function(x,...){
    attr(x,"class")<-NULL
    attr(x,"layout")<-NULL
    print(x)
}
