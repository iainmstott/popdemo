# helpers for project.compadreXXX

# switch for aPreChecks and a few others
.isDeterministic <- function(aMat) {
  (is.list(aMat) & length(aMat) == 1) |
    is.matrix(aMat)
}

# Blanket function to cover both stochastic and deterministic
.aPreChecks <- function(aMat, errors) {
  
  if(.isDeterministic(aMat)) {
    errors <- c(errors, .checkDeterministicA(aMat, errors))
    
  } else {
    
    errors <- c(errors, .checkStochasticA(aMat, errors))
  }
  
  return(errors)
}

.isSquare <- function(aMat) {
  dim(aMat)[1] == dim(aMat)[2]
}

# check a single matrix for deterministic iterations
.checkDeterministicA <- function(aMat, errors) {
  
  aMat <- aMat[[1]]  
  if(!.isSquare(aMat)) {
    errors <- c(errors, 'A must be a square matrix')
  }
  
}

# converts A slot of CompadreMSM to an array so that Iain's original algorithm
# works as is.
.a2Array <- function(aList) {
  
  if(.isDeterministic(aList)) {
    if(is.list(aList)) A <- aList[[1]]
    if(is.matrix(aList)) A <- aList
    M1 <- A
    dim(A) <- c(dim(A), 1)
    dimnames(A)[[1]] <- dimnames(M1)[[1]]
    dimnames(A)[[2]] <- dimnames(M1)[[2]]
    dimnames(A)[[3]] <- NULL  
  } else {
    
    numA <- length(aList)
    dimA <- dim(aList[[1]])[1]
    
    A <- numeric(dimA * dimA * numA)
    dim(A) <- c(dimA, dimA, numA)
    
    # generate array of As
    for(i in seq_len(numA)){
      A[,,i] <- aList[[i]]
    }
    
    dimnames(A)[[1]] <- dimnames(aList[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(aList[[1]])[[2]]
    dimnames(A)[[3]] <- names(aList)
  }
  
  return(A)
}



# Check set of matrices for stochastic iterations
.checkStochasticA <- function(aMat, errors) {
  
  # stop if A isn't a matrix or list of matrices (or array of matrices)
  if(!any((is.list(aMat) & all(vapply(aMat, is.matrix, logical(1)))), 
          (is.array(aMat) & length(dim(aMat)) == 3)) ){
    c(errors, "A must be a matrix or list of matrices")
  }
  
  # test all equal dimension
  allDim <- vapply(aMat, dim, numeric(2))
  
  if(!diff(range(allDim)) == 0) {
    c(errors, "all matrices in A must be square and have the same dimension as each other")
  }
  # test squareness
  squares <- vapply(aMat, .isSquare, logical(1))
  if(any(!squares)) {
    squareInd <- which(!squares)
    errors <- c(errors, 
                paste(c('One or more matrices are not square (',
                        paste(squareInd, 
                              collapes = ', '),
                        ')'),
                      collapse = ""))
  }
  
  return(errors)
}

# only used if(PREcheck)
.checkRedAndPrim <- function(A) {
  # optional primitivity/reducibility checks
  
  if(length(A) > 1){
    reds <- vapply(A, isIrreducible, logical(1))
    imps <- vapply(A, isPrimitive, logical(1))
    
    if(any(!reds)) {
      redInd <- which(!reds)
      warning(paste(c("One or more matrices are reducible (", 
                      paste(redInd, collapse = ", "),")"),
                    sep = ""))
    }
    
    if(any(!imps)) {
      impInd <- which(!imps)
      warning(paste(c("One or more matrices are reducible (", 
                      paste(impInd, collapse = ", "),")"),
                    sep = ""))
    }
  }
  
  if(length(A) == 1) {
    A <- A[[1]]
    if (!isIrreducible(A)) {
      warning("Matrix is reducible")
    } else {
      if (!isPrimitive(A)) {
        warning("Matrix is imprimitive")
      }
    }
  }
  
  invisible(NULL)
}

# checks aStart for most problems. Generates MCstart if none are found
.checkAstart <- function(Astart, stochSeq, matrixNames, nMat) {
  
  MCstart <- NULL
  
  
  if(nMat == 1) {
    if(stochSeq != "unif") warning("Ignoring stochSeq: only 1 matrix in A")
  }
  
  if(nMat > 1) {
    
    # check types
    if(!any(is.integer(Astart),
            is.character(Astart),
            is.null(Astart))) stop("Astart should be numeric, character or NULL")
    
    # check all numerics are positive
    if(is.integer(Astart)) {
      if(length(Astart) > 1) stop("Astart should be length 1")
      if(Astart < 1) stop("Astart must be greater than 0")
      if(Astart > nMat) stop("Astart cannot be greater than the number of matrices in A")
      if(Astart%%1 != 0) stop("Astart must be an integer")
      MCstart <- Astart
    }
    if(is.character(Astart)) {
      if(length(Astart) > 1) stop("Astart should be length 1")
      if(!(Astart %in% matrixNames)) stop("Astart is not found in names of matrices in A")
      MCstart <- match(Astart, matrixNames)
    }
  }
  
  # if it's null, pick a random one
  if(is.null(Astart)) {
    matInds <- rmultinom(1, 1, rep(1/nMat, nMat))
    MCstart <- which(matInds == 1)
  }
  
  return(MCstart)
}

# checks stochSeq. Does not return anything if there no problems because
# stochSeq is then passed as is to initMc by checkAstart(). 
# 
# This is a mess of if()s, but I'm not sure how to avoid this.

.checkstochSeq <- function(stochSeq, matrixNames, matDim) {
  
  # general checks
  if(!any(stochSeq[1] == "unif", 
          is.matrix(stochSeq), 
          is.integer(stochSeq) & is.null(dim(stochSeq)),
          is.character(stochSeq) & is.null(dim(stochSeq)))){
    stop('stochSeq must take "unif", a numeric matrix, a numeric vector, or a character vector',
         call. = FALSE)
  }
  
  # if it's a transition matrix
  if(is.matrix(stochSeq)) {
    if(dim(stochSeq)[1] != length(matrixNames)) {
      stop('Dimension of stochSeq must equal the number of matrices in A',
           call. = FALSE)
    }
    if(dim(stochSeq)[1] != dim(stochSeq)[2]) {
      stop('stochSeq is not square', call. = FALSE)
    }
    if(!all(colSums(stochSeq) == 1)) {
      stop("Columns of stochSeq do not sum to 1")
    }
  }
  
  # if it's an integer sequence 
  if(is.integer(stochSeq)) {
    if(min(stochSeq) < 1) stop('All entries in stochSeq must be positive integers')
    if(max(stochSeq) > length(matrixNames)) {
      stop('One or more entries in stochSeq are greater than the number of matrices',
           'in A')
    }
  }
  
  # if it's a character sequence w/ matrix names
  if(is.character(stochSeq) &
     length(stochSeq) > 1 &
     !all(stochSeq %in% matrixNames)) {
    stop("Names of stochSeq aren't all found in A")  
  }
  
  invisible(NULL)
}

# takes the user-specified input and generates a sequence of integers to index
# A by. This index is used to select the A matrix during iteration
.initMc <- function(stochSeq, time, MCstart, matrixNames, nMat) {
  # browser()
  
  # if uniform prob, create a uniform transition matrix
  # if user provided markov matrix, just use that to generate markov chains
  if(stochSeq == "unif" || is.matrix(stochSeq)) {
    
    # user supplied matrix
    if(is.matrix(stochSeq)) {
      MCtm <- stochSeq 
      
    } else {
      # Unif gets uniform probability matrix
      MCtm <- matrix(rep(1/nMat, nMat^2), nMat, nMat)
    }
    
    # convert matrix to index
    MC <- .rmc(MCtm, time, MCstart)
    
  } else if(is.integer(stochSeq) & is.null(dim(stochSeq))) {
    MC <- stochSeq
  } else if(is.character(stochSeq)) {
    MC <- match(stochSeq, matrixNames)
  }
  
  if(!is.null(MCstart)) {
    MC[1] <- MCstart
  }
  
  names(MC) <- matrixNames[MC]
  
  return(MC)
}

#' @noRd
# takes user-sepcified input and generates an initial population vector or vectors
# for usage in iteration. Always 3 dimensions
# a given column is for each stage class, rows are time,
# and the third dimension is for each set of initial vectors

.initPopVec <- function(A, 
                        vector,
                        time, 
                        matDim, 
                        stageNames,
                        draws, 
                        alpha.draws,
                        standard.vec) {
  # browser()
  # new character variable for switching between initial vectors
  if(is.numeric(vector)) {
    vectorSwitch <- 'user'
  } else {
    vectorSwitch <- vector
  }
  
  popVec <- switch(vectorSwitch,
                   'diri' = .initDiriVector(matDim = matDim,
                                            time = time,
                                            stageNames = stageNames,
                                            draws = draws,
                                            alpha.draws = alpha.draws),
                   'n' = .initBiasVector(matDim = matDim,
                                         time = time,
                                         stageNames = stageNames),
                   'user' = .initUserVector(vector = vector,
                                            time = time,
                                            matDim = matDim, 
                                            stageNames = stageNames,
                                            standard.vec = standard.vec),
                   # not yet implemented, waiting to see what
                   # the database version of these looks like
                   # 'db' = .initDbVector(vector, 
                   #                      matDim, 
                   #                     stageNames)
  )
  
  # can be any combinatin of c(deterministic, stochastic) AND
  # c(single, multiple, diri, or bias)
  projType <- ifelse(dim(A)[[3]] > 1, 'stochastic', 'deterministic')
  if(is.numeric(vector)) {
    if((is.array(popVec) & dim(popVec)[3] == 1)) vecType <- 'single'
    if(dim(popVec)[3] > 1) vecType <- 'multiple'
    
  } else {
    vecType <- vector
  }
  
  out <- list(popVec = popVec,
              projType = projType,
              vecType = vecType)
  
  return(out)
}

#' @noRd
.initBiasVector <- function(matDim, time, stageNames) {
  
  
  VecBias <- numeric((time + 1) * matDim * matDim)
  dim(VecBias) <- c(time + 1, matDim, matDim)
  dimnames(VecBias)[[2]] <- stageNames
  dimnames(VecBias)[[3]] <- paste("bias", stageNames, sep = "")
  
  
  # initial vectors come from identity matrix
  I <- diag(matDim)
  VecBias[1, ,] <- I
  
  return(VecBias)
}

#' @noRd
.initDiriVector <- function(matDim, time, stageNames, draws, alpha.draws) {
  
  # dirichlet parameters
  if(alpha.draws[1] == "unif") alpha.draws <- rep(1, matDim)
  # 
  if(length(alpha.draws) != matDim) {
    stop("length of alpha.draws must be equal to matrix dimension")
  }
  
  # array of Dirichlet vectors
  # a given column is for each stage class, rows are time,
  # and the third dimension is for each set of possible initial conditions
  VecDraws <- numeric((time + 1) * matDim * draws)
  dim(VecDraws) <- c(time + 1, matDim, draws)
  dimnames(VecDraws)[[2]] <- stageNames
  dimnames(VecDraws)[[3]] <- paste("draw", as.character(1:draws), sep = "")
  
  # initial stage vectors come from dirichlet
  VecDraws[1, , ] <- t(MCMCpack::rdirichlet(draws, alpha.draws))
  
  return(VecDraws)
}

.initUserVector <- function(vector, time, matDim, stageNames, standard.vec) {
  # single vector
  if(length(vector) == matDim) {
    
    # standardisation
    if (standard.vec) vector <- .standardVec(vector)
    # empty objects
    Vec <- array(0, c(time + 1, matDim, 1))
    dimnames(Vec)[[2]] <- stageNames
    
    Vec[1, , ] <- vector
    
    dimnames(Vec)[[3]] <- paste("V", 1:dim(Vec)[3], sep = "")
    
  } else { # multiple vectors
    
    if(length(vector) %% matDim != 0) stop('Vector has wrong dimensions')
    
    nvec <- dim(vector)[2] 
    
    vectornames <- dimnames(vector)[[2]]
    if(is.null(vectornames)) {
      vectornames <- paste("V", as.character(1:nvec), sep="")
    }
    
    # standardisations
    if (standard.vec) vector <- apply(vector, 2, .standardVec)
    
    Vec <- numeric((time + 1) * matDim * nvec)
    dim(Vec) <- c(time + 1, matDim, nvec)
    dimnames(Vec)[[2]] <- stageNames
    dimnames(Vec)[[3]] <- vectornames
    Vec[1, , ] <- vector
  }
  
  return(Vec)
}

.standardVec <- function(vec) {
  
  vec/sum(vec)
}

#' @noRd
.initPopSum <- function(nVec, time, vector) {
  popSize <- numeric((time + 1) * nVec)
  dim(popSize) <- c(time + 1, nVec)
  
  # The vectors are standardized (or not) when generated, so we don't need
  # to worry about accounting for those here
  if(dim(vector)[[3]] == 1) {
    popSize[1, ] <- rowSums(vector)[1]
    
    # Copy over dimnames for plotting
    dimnames(popSize)[[2]] <- list(dimnames(vector)[[3]])
    
  } else {
    
    popSize[1, ] <- apply(vector, 3, rowSums)[1, ]
    
    dimnames(popSize)[[2]] <- dimnames(vector)[[3]]
  }
  
  return(popSize)
}

# does the iterating
.iterateMsm <- function(A, time, 
                        vecs, popSize,
                        MC, 
                        vecType, projType,
                        return.vec) {
  vecOut <- vecs
  popOut <- popSize
  nReps <- dim(vecOut)[[3]]
  
  # iterates over the different vectors. If single vector, then break after i = 1
  for(i in 1:nReps) {
    
    for(j in seq_len(time)){
      vecOut[j + 1, , i] <- A[, , MC[j]] %*% vecOut[j, , i]
      popOut[j + 1, i] <- sum(vecOut[j + 1, , i])
    }
    
    
    if(projType == 'deterministic') {
      
      bounds <- t(apply(popOut, 1, range))
      
    }  else {
      # no bounds for stochastic projections
      bounds <- matrix(c(NA_real_, NA_real_), nrow = 1)
    }
  }
  
  out <- .updateProjOutput(vecOut,
                           popOut,
                           bounds, 
                           A, 
                           time,
                           MC, 
                           vecType, 
                           projType, 
                           return.vec)
  
  return(out)
}

#' Standardize a single or set of matrices
#' 
#' Standardizes matrices in an array by their dominant eigenvalues. This filters
#' out the asymptotic components of each matrix, leaving only the transient
#' growth components.
#' 
#' @param A A single square matrix or array of square matrices stacked as a
#' 3-dimensional array
#' 
#' @return An array of the same dimension as the input (\code{A}). They are
#' the input matrices that have been divided by their dominant eigenvalues.
#' 
#' @references Add something here
#' 
#' @export

standardA <- function(A) {
  matDim <- dim(A)[1]
  nMat <- dim(A)[3]
  # placeholder
  Ms <- A
  # get lambdas
  lambda <- apply(Ms, 3, eigs, what = "lambda")
  # initialize output
  A <- numeric(matDim * matDim * nMat)
  dim(A) <- c(matDim, matDim, nMat)
  
  # loop to get standardized As
  for(i in seq_len(nMat)){
    A[,,i] <- Ms[,,i]/lambda[i]
  }
  
  return(A)
}

.a2List <- function(A) {
  out <- list()
  nMat <- dim(A)[[3]]
  
  for(i in seq_len(nMat)) {
    out[[i]] <- A[ , , i]
  }
  
  names(out) <- dimnames(A)[[3]]
  
  return(out)
  
}


#' @noRd
# updates Projection components at the end of each iteration
.updateProjOutput <- function(vecOut,
                              popOut,
                              bounds, 
                              A,
                              time,
                              MC,
                              vecType,
                              projType,
                              return.vec) {
  
  out <- Projection()
  if(is.null(dim(popOut))) dim(popOut) <- time + 1
  out@.Data <- popOut
  if(return.vec) out@vec <- vecOut
  if(is.list(A)) {
    out@mat <- A
  } else {
    out@mat <- .a2List(A)
  }
  out@stochSeq <- as.integer(MC)
  out@projtype <- projType
  out@bounds <- bounds
  out@vectype <- vecType
  
  return(out) 
}


.errorConstructor <- function(errors) {
  n <- seq(1, length(errors), 1)
  
  paste(
    c(
      'The following errors were detected:\n',
      paste(n, 
            '. ',
            errors,
            '\n'),
      sep = ""
    ),
    sep = ""
  )
}
