vec2var.jomoni <- function (z, r = 1){
  if (!(class(z) == "ca.jomoni")) {
    stop("\nPlease, provide object of class 'ca.jomoni' as 'z'.\n")
  }
  r <- as.integer(r)
  if (!({
    1 <= r
  } && {
    r < ncol(z@x)
  })) {
    stop(paste("\nThe cointegration rank 'r' must be in the interval [2:", 
               ncol(z@x) - 1, "].\n", sep = ""))
  }
  etc <- z@ZK %*% z@V[, 1:r]
  colnames(etc) <- paste("etc", 1:r, sep = "")
  coeffs <- coef(lm(z@Z0 ~ -1 + as.numeric(etc) + z@Z1))
  rownames(coeffs) <- c(colnames(etc), colnames(z@Z1))
  PI <- z@W[, 1:r] %*% t(z@V[, 1:r])
  # Number of columns in the break.matrix
  n.breaks <- ncol(z@break.matrix)
  if (is.null(n.breaks)){
    n.breaks <- 0
  } else {
    n.breaks <- as.integer(ncol(z@break.matrix))
  }
  
  # Additional if conditions to account for all cases possible with ca.jomoni
  # Take into account the number of cols, which are part of the break.matrix
  # -  no constant  --nc
  # -  restricted constant, i.e. constant in the cointegration vector --rc
  # -  unrestricted constant, i.e. constant in the VECM equation, not in the cointegration vector --uc
  # -  constant + restricted trend, i.e. constant in the VECM equation and restricted trend in the cointegration vector --crt
  # -  constant + unrestricted trend, i.e. constant and trend in the VECM equation, not in the cointegration vector --ct
  if (z@ecdet == "nc") {
    if(n.breaks == 0){
      detcoeffs <- matrix( nrow = 0, ncol = ncol(z@x), 
                           byrow = TRUE)
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI
      rhs <- z@Z1
      } else {
      detcoeffs <- matrix(PI[, (z@P + 1):(z@P + n.breaks)], nrow = n.breaks, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, (z@P + 1): (z@P + n.breaks)])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -((z@P + 1): (z@P + n.breaks))]
      rhs <- z@Z1 
      # Add the columns of the break.matrix
      # Add the additional constant or trend
      rhs <- cbind(z@break.matrix[-(1:z@lag),], z@Z1)
      colnames(rhs) <- c(colnames(z@break.matrix), colnames(z@Z1))
    }
    
  } else if (z@ecdet == "rc") {  
    if(n.breaks == 0){
      detcoeffs <- matrix(PI[, z@P + 1], nrow = 1, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- "constant"
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -(z@P + 1)]
      rhs <- cbind(1, z@Z1)
      colnames(rhs) <- c("constant", colnames(z@Z1))
      
    } else {
      detcoeffs <- matrix(PI[, (z@P + 1) : (z@P + n.breaks + 1) ], nrow = n.breaks + 1, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, (z@P + 1) : (z@P +  n.breaks + 1)])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -((z@P + 1) : (z@P +  n.breaks + 1))]
      rhs <-  z@Z1
      rhs <- cbind(1, z@break.matrix[-(1:z@lag),], z@Z1)
      colnames(rhs) <- c("constant", colnames(z@break.matrix), colnames(z@Z1))
    }
    
  } else if (z@ecdet == "uc") { 
    if(n.breaks == 0){
      detcoeffs <- matrix( nrow = 0, ncol = ncol(z@x), 
                           byrow = TRUE)
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI
      rhs <- z@Z1
    } else {
      detcoeffs <- matrix(PI[, (z@P + 1) : (z@P +  n.breaks)], nrow = n.breaks, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, (z@P + 1) : (z@P +  n.breaks)])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -((z@P + 1) : (z@P +  n.breaks))]
      rhs <- cbind(z@break.matrix[-(1:z@lag),], z@Z1)
      colnames(rhs) <- c(colnames(z@break.matrix), colnames(z@Z1))
      
    }
    
  } else if (z@ecdet == "crt") {
    if(n.breaks == 0){
      detcoeffs <- matrix(PI[, z@P + 1], nrow = 1, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, z@P + 1])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -(z@P + 1)]
      
      rhs <- cbind( z@ZK[, z@P + 1], z@Z1)
      colnames(rhs) <- c(colnames(z@ZK)[z@P+1], colnames(z@Z1))
      
    } else {
      detcoeffs <- matrix(PI[, (z@P + 1) : (z@P + n.breaks + 1) ], nrow = n.breaks + 1, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, (z@P + 1) : (z@P +  n.breaks + 1)])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -((z@P + 1) : (z@P +  n.breaks + 1))]
      rhs <- cbind( z@ZK[, z@P + 1], z@break.matrix[-(1:z@lag),], z@Z1)
      colnames(rhs) <- c(colnames(z@ZK)[z@P+1], colnames(z@break.matrix), colnames(z@Z1))
      
    }
  } else if (z@ecdet == "ct"){
    if(n.breaks == 0){
      detcoeffs <- matrix( nrow = 0, ncol = ncol(z@x), 
                           byrow = TRUE)
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI
      rhs <- z@Z1
    } else {
      detcoeffs <- matrix(PI[, (z@P + 1) : (z@P +  n.breaks)], nrow = n.breaks, ncol = ncol(z@x), 
                          byrow = TRUE)
      rownames(detcoeffs) <- colnames(PI[, (z@P + 1) : (z@P +  n.breaks)])
      colnames(detcoeffs) <- colnames(z@x)
      PI <- PI[, -((z@P + 1) : (z@P +  n.breaks))]
      rhs <- cbind(z@break.matrix[-(1:z@lag),], z@Z1)
      colnames(rhs) <- c(colnames(z@break.matrix), colnames(z@Z1))
    }
  }
  
  if (!(is.null(eval(z@season)))) {
    seas <- eval(z@season) - 1
    season <- paste("sd", 1:seas, sep = "")
    detcoeffs <- rbind(detcoeffs, coeffs[season, ])
  }
  if (!(is.null(eval(z@dumvar)))) {
    dumnames <- colnames(z@dumvar)
    tmp <- rownames(detcoeffs)
    detcoeffs <- rbind(detcoeffs, coeffs[dumnames, ])
    rownames(detcoeffs) <- c(tmp, dumnames)
  }
  detcoeffs <- t(detcoeffs)
  Gamma <- t(coeffs[-which(rownames(coeffs) %in% c(colnames(detcoeffs), 
                                                   colnames(etc))), ])
  rownames(Gamma) <- colnames(z@x)
  A <- list()
  if (identical(z@spec, "transitory")) {
    if (identical(z@lag, as.integer(2))) {
      A$A1 <- Gamma + PI + diag(z@P)
      rownames(A$A1) <- colnames(z@x)
      colnames(A$A1) <- paste(colnames(z@x), ".l1", sep = "")
      A$A2 <- -1 * Gamma
      rownames(A$A2) <- colnames(z@x)
      colnames(A$A2) <- paste(colnames(z@x), ".l2", sep = "")
    }
    else if (z@lag > 2) {
      idx.end <- seq(from = z@P, by = z@P, length.out = z@lag - 
                       1)
      idx.start <- idx.end - z@P + 1
      A[[1]] <- Gamma[, idx.start[1]:idx.end[1]] + PI + 
        diag(z@P)
      rownames(A[[1]]) <- colnames(z@x)
      colnames(A[[1]]) <- paste(colnames(z@x), ".l1", sep = "")
      for (i in 2:(z@lag - 1)) {
        A[[i]] <- Gamma[, idx.start[i]:idx.end[i]] - 
          Gamma[, idx.start[i - 1]:idx.end[i - 1]]
        rownames(A[[i]]) <- colnames(z@x)
        colnames(A[[i]]) <- paste(colnames(z@x), ".l", 
                                  i, sep = "")
      }
      A[[z@lag]] <- -1 * Gamma[, tail(idx.start, 1):tail(idx.end, 
                                                         1)]
      rownames(A[[z@lag]]) <- colnames(z@x)
      colnames(A[[z@lag]]) <- paste(colnames(z@x), ".l", 
                                    z@lag, sep = "")
      names(A) <- paste("A", 1:z@lag, sep = "")
    }
  }
  if (identical(z@spec, "longrun")) {
    if (identical(z@lag, as.integer(2))) {
      A$A1 <- Gamma + diag(z@P)
      rownames(A$A1) <- colnames(z@x)
      colnames(A$A1) <- paste(colnames(z@x), ".l1", sep = "")
      A$A2 <- PI + diag(z@P) - A$A1
      rownames(A$A2) <- colnames(z@x)
      colnames(A$A2) <- paste(colnames(z@x), ".l2", sep = "")
    }
    else if (z@lag > 2) {
      idx.end <- seq(from = z@P, by = z@P, length.out = z@lag - 
                       1)
      idx.start <- idx.end - z@P + 1
      A[[1]] <- Gamma[, idx.start[1]:idx.end[1]] + diag(z@P)
      rownames(A[[1]]) <- colnames(z@x)
      colnames(A[[1]]) <- paste(colnames(z@x), ".l1", sep = "")
      for (i in 2:(z@lag - 1)) {
        A[[i]] <- Gamma[, idx.start[i]:idx.end[i]] - 
          Gamma[, idx.start[i - 1]:idx.end[i - 1]]
        rownames(A[[i]]) <- colnames(z@x)
        colnames(A[[i]]) <- paste(colnames(z@x), ".l", 
                                  i, sep = "")
      }
      A[[z@lag]] <- PI - Gamma[, tail(idx.start, 1):tail(idx.end, 
                                                         1)]
      rownames(A[[z@lag]]) <- colnames(z@x)
      colnames(A[[z@lag]]) <- paste(colnames(z@x), ".l", 
                                    z@lag, sep = "")
      names(A) <- paste("A", 1:z@lag, sep = "")
    }
  }
  datamat <- embed(z@x, dimension = z@lag + 1)
  datamat <- cbind(datamat[, 1:ncol(z@x)], rhs[, colnames(detcoeffs)], 
                   datamat[, -c(1:ncol(z@x))])
  temp1 <- NULL
  for (i in 1:z@lag) {
    temp <- paste(colnames(z@x), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  colnames(datamat) <- c(colnames(z@x), colnames(detcoeffs), 
                         temp1)
  resids <- datamat[, colnames(z@x)] - datamat[, colnames(detcoeffs)] %*% 
    t(detcoeffs)
  for (i in 1:z@lag) {
    resids <- resids - datamat[, colnames(A[[i]])] %*% t(A[[i]])
  }
  colnames(resids) <- paste("resids of", colnames(z@x))
  result <- list(deterministic = detcoeffs, A = A, p = z@lag, 
                 K = ncol(z@x), y = z@x, obs = nrow(z@Z0), totobs = nrow(z@x), 
                 call = match.call(), vecm = z, datamat = datamat, resid = resids, 
                 r = r)
  class(result) <- "vec2var"
  return(result)
}