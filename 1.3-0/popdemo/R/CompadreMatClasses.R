#' @rdname Stochastic-Matrix-Classes
#' @title Stochastic Matrix Classes 
#' 
#' @description Stochastic matrix classes in \code{Rcompadre} and \code{popdemo} can be broken
#' down into two type: stochastic by matrix selection and stochastic by element
#' selection. Within the latter category, there are three different types supported:
#' density dependent, environmentally dependent, and their combination. 
#' 
#' 

CompadreDDM <- setClass('CompadreDDM', contains = "list",
                        slots = c(data_list = 'list',
                                  mat_exprs = 'list'))

ValidCompadreDDM <- function(object) {
  

  if(!'mat_expr' %in% names(object@mat_exprs)) {
    stop('User-defined "CompadreDDM" objects must contain a "mat_expr" expression.\n',
         'See ?make_mat_exprs for more details.')
  }
  
  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreDDM", ValidCompadreDDM)

CompadreEPM <- setClass('CompadreEPM', contains = "list",
                        slots = c(data_list = 'list',
                                  mat_exprs = 'list',
                                  env_exprs = 'list'))


ValidCompadreEPM <- function(object) {

  if(!'mat_expr' %in% names(object@mat_exprs)) {
    stop('User-defined "CompadreEPM" objects must contain a "mat_expr" expression.\n',
         'See ?make_mat_exprs for more details.')
  }
  
  # need to add ability to check for potential negative matrix values.
  TRUE
}

setValidity("CompadreEPM", ValidCompadreEPM)

#' @rdname Stochastic-Matrix-Classes
#' @export
setGeneric("mat_expr", function(object) {
  standardGeneric('mat_expr')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_expr', signature = c("CompadreDDM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@mat_exprs$mat_expr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_expr', signature = c("CompadreEPM"),
          function(object) {
            temp <- rlang::quo_get_expr(object@mat_exprs$mat_expr)
            # out <- expr_to_char_mat(temp) # need to write this
            return(temp)
          }
)

#' @rdname Stochastic-Matrix-Classes
#' @export
setGeneric('mat_funs', function(object) {
  standardGeneric('mat_funs')
})

#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_funs', signature = c("CompadreDDM"),
          function(object) {
            ind <- which(names(object@mat_exprs) == 'mat_expr')
            funs <- object@mat_exprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


#' @rdname Stochastic-Matrix-Classes
#' @export
setMethod('mat_funs', signature = c("CompadreEPM"),
          function(object) {
            ind <- which(names(object@mat_exprs) == 'mat_expr')
            funs <- object@mat_exprs[-ind]
            out <- lapply(funs, rlang::quo_get_expr)
            return(out)
          }
)


