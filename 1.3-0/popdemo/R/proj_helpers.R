# helpers for project.compadreXXX

check_a_start <- function(a_start) {
  
  
}
  
  
check_a_seq <- function(a_seq, a_seq_type) {
  
  
}

a_pre_checks <- function(a_mat) {
  if (!isIrreducible(a_mat[,,1])) {
    warning("Matrix is reducible")
  } else {
    if (!isPrimitive(A[,,1])) {
      warning("Matrix is imprimitive")
    }
  }
  
}

check_a_mat <- function(a_list) {
  
  
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
