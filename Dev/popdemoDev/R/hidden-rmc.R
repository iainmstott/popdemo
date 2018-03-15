#' @importFrom stats runif rmultinom
.rmc <- function(tm, cl, s1 = NULL){
    tm_dim <- dim(tm)
    if(length(tm_dim) != 2 | tm_dim[1] != tm_dim[2]){
        stop("Markov transition matrix must be square")
    }
    if(!all(colSums(tm)==1)){
        stop("all column sums of Markov transition matrix must equal 1")
    }
    tm_cum <- apply(tm, 2, cumsum)
    if(!is.null(s1)){
        if(any(!(s1%in%(1:tm_dim[1])), s1%%1 != 0)){
            stop("starting state (s1) must be integer with 0 < s1 <= S, \n where S is Markov transition matrix dimension")
        }
    }
    states <- numeric(cl)
    ifelse(is.null(s1),
           states[1] <- which(stats::rmultinom(1, 1, rep(1, dim(tm)[1])) == 1),
           states[1] <- s1
    )
    rd <- stats::runif(cl - 1)
    for(i in 2:cl) states[i] <- sum(rd[i-1] > tm_cum[, states[i-1]]) + 1
    return(states)
}
