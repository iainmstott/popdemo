# helpers for project.compadreXXX

is_deterministic <- function(A) {
  (is.list(A) & length(A) == 1) |
    is.matrix(A)
}


check_a_start <- function(a_start, errors) {
  
  
}
  
  
check_a_seq <- function(a_seq, a_seq_type, errors) {
  
  
}

a_pre_checks <- function(a_mat, errors, PREcheck) {
  
  if(is_deterministic(A)) {
    errors <- c(errors, check_deterministic_A(a_mat, errors, PREcheck))
    
  } else {
  
    errors <- c(errors, check_stochastic_A(a_mat, errors, PREcheck))
  }
  
  return(errors)
}

is_square <- function(a_mat) {
  dim(a_mat)[1] == dim(a_mat)[2]
}

check_deterministic_A <- function(a_mat, errors, PREcheck) {
  
  if(!is_square(a_mat)) {
    errors <- c(errors, 'A must be a square matrix')
  }
  
  if(PREcheck){
    if (!isIrreducible(a_mat)) {
      warning("Matrix is reducible")
    } else {
      if (!isPrimitive(a_mat)) {
        warning("Matrix is imprimitive")
      }
    }
  }
  
}

a_to_array <- function(a_list) {
  
  if(is_deterministic(a_list)) {
    if(is.list(a_list)) A <- a_list[[1]]
    if(is.matrix(a_list)) A <- a_list
    M1 <- A
    dim(A) <- c(dim(A), 1)
    dimnames(A)[[1]] <- dimnames(M1)[[1]]
    dimnames(A)[[2]] <- dimnames(M1)[[2]]
    dimnames(A)[[3]] <- NULL  
  } else {
    
    numA <- length(a_list)
    dimA <- dim(a_list[[1]])[1]
    
    A <- numeric(dimA * dimA * numA)
    dim(A) <- c(dimA, dimA, numA)
    
    # generate array of As
    for(i in seq_len(numA)){
      A[,,i] <- a_list[[i]]
    }
    
    dimnames(A)[[1]] <- dimnames(a_list[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(a_list[[1]])[[2]]
    dimnames(A)[[3]] <- names(a_list)
  }
}

check_stochastic_A <- function(a_mat, errors, PREcheck) {
  
  # stop if A isn't a matrix or list of matrices (or array of matrices)
  if(!any((is.list(a_mat) & all(vapply(a_mat, is.matrix, logical(1)))), 
          (is.array(a_mat) & length(dim(a_mat)) == 3)) ){
    c(errors, "A must be a matrix or list of matrices")
  }
  
  # test all equal dimension
  all_dim <- vapply(a_mat, dim, numeric(2))
  if(!diff(range(all_dim)) == 0) {
    c(errors, "all matrices in A must be square and have the same dimension as each other")
  }
  # test squareness
  squares <- vapply(a_mat, is_square, logical(1))
  if(!any(squares)) {
    errors <- c(error, 'All matrices in A must be square')
  }
  
  # optional primitivity/reducibility checks
  if(PREcheck){
    reds <- vapply(a_mat, isIrreducible, logical(1))
    imps <- vapply(a_mat, isPrimitive, logical(1))
    
    if(any(!reds)) {
      red_ind <- which(!reds)
      warning(paste(c("One or more matrices are reducible (", 
                      paste(red_ind, collapse = ", "),")"),
                    sep = ""))
    }
    
    if(any(!imps)) {
      imp_ind <- which(!imps)
      warning(paste(c("One or more matrices are reducible (", 
                      paste(imp_ind, collapse = ", "),")"),
                    sep = ""))
    }
    
  }
  
  return(errors)
}
init_pop_vec <- function(vector, ...) {
  
  # new character variable for switching between initial vectors
  if(is.numeric(vector)) {
    vector_switch <- 'user'
  } else {
    vector_switch <- vector
  }
  
  switch(vector_switch,
         'diri' = init_diri_vector(n_classes, alpha.draws),
         'n' = init_bias_vector(n_classes),
         'user' = init_user_vector(vector),
         'db' = init_db_vector(vector, n_classes, alpha.draws))
}

init_a_seq <- function(a_seq) {
  init_unif()
  
  init_mc()
  
  init_numeric()
  
  init_char()
}
  
iterate_msm <- function(...) {
  
}


initialize_proj_outputs <- function(n_generations,
                                    vector,
                                    init_p_vec) {
  
  output <- Projection()
  output@.Data <- sum(init_p_vec)
  
  # Population vectors
  output@vec <- array(
    t(init_p_vec),
    dim = c(1, 
            ifelse(is.numeric(vector) & !is.matrix(vector), 
                   length(init_p_vec),
                   dim(init_p_vec)[1]),
            ifelse(is.numeric(vector) & !is.matrix(vector), 2, 3)
    )
  )
  
  
  
  if('growth_rates' %in% target_output) {
    output$growth_rates <- data.frame(n_tot = c(sum(init_p_vec),
                                                rep(NA_real_, n_generations - 1)),
                                      lambda = rep(NA_real_, n_generations))
  }
  
  return(output)
  
}

update_proj_outputs <- function(object, to_insert, iteration) {
  
  object@.Data[iteration] <- sum(to_insert)
  object@vec[iteration] <- t(to_insert)
  
}
