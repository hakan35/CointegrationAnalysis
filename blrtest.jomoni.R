
blrtest.jomoni <- function (z, H, r){
  if (!(class(z) == "ca.jomoni")) {
    stop("\nPlease, provide object of class 'ca.jomoni' as 'z'.\n")
  }
  if (r >= z@P || r < 1) {
    stop("\nCount of cointegrating relationships is out of allowable range.\n")
  }
  # Take into account the different possible options of the ca.jomoni function
  # c("nc", "rc", "uc", "crt", "ct")
  # Include the possibility to define one of the five cases for the VECM
  # -  no constant  --nc
  # -  restricted constant, i.e. constant in the cointegration vector --rc
  # -  unrestricted constant, i.e. constant in the VECM equation, not in the cointegration vector --uc
  # -  constant + restricted trend, i.e. constant in the VECM equation and restricted trend in the cointegration vector --crt
  # -  constant + unrestricted trend, i.e. constant and trend in the VECM equation, not in the cointegration vector --ct
  if (z@ecdet == "nc" || z@ecdet == "uc" || z@ecdet == "ct") {
    P <- z@P + ncol(z@break.matrix)
  }
  else {
    P <- z@P + 1 + ncol(z@break.matrix)
  }
  r <- as.integer(r)
  H <- as.matrix(H)
  if (!(nrow(H) == P)) {
    stop("\nRow number of 'H' is unequal to VAR order.\n")
  }
  type <- "Estimation and testing under linear restrictions on beta"
  N <- nrow(z@Z0)
  M00 <- crossprod(z@Z0)/N
  M11 <- crossprod(z@Z1)/N
  MKK <- crossprod(z@ZK)/N
  M01 <- crossprod(z@Z0, z@Z1)/N
  M0K <- crossprod(z@Z0, z@ZK)/N
  MK0 <- crossprod(z@ZK, z@Z0)/N
  M10 <- crossprod(z@Z1, z@Z0)/N
  M1K <- crossprod(z@Z1, z@ZK)/N
  MK1 <- crossprod(z@ZK, z@Z1)/N
  M11inv <- solve(M11)
  S00 <- M00 - M01 %*% M11inv %*% M10
  S0K <- M0K - M01 %*% M11inv %*% M1K
  SK0 <- MK0 - MK1 %*% M11inv %*% M10
  SKK <- MKK - MK1 %*% M11inv %*% M1K
  Ctemp <- chol(t(H) %*% SKK %*% H, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  S00inv <- solve(S00)
  valeigen <- eigen(Cinv %*% t(H) %*% SK0 %*% S00inv %*% S0K %*% 
                      H %*% t(Cinv))
  e <- valeigen$vector
  V <- H %*% t(Cinv) %*% e
  Vorg <- V
  idx <- ncol(V)
  V <- sapply(1:idx, function(x) V[, x]/V[1, x])
  W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)
  PI <- W %*% t(V)
  DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% 
    t(V) %*% SK0
  GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv
  lambda.res <- valeigen$values
  lambda <- z@lambda
  teststat <- N * sum(log((1 - lambda.res[1:r])/(1 - lambda[1:r])))
  df <- r * (P - ncol(H))
  pval <- c(1 - pchisq(teststat, df), df)
  new("cajo.test", Z0 = z@Z0, Z1 = z@Z1, ZK = z@ZK, ecdet = z@ecdet, 
      H = H, A = NULL, B = NULL, type = type, teststat = teststat, 
      pval = pval, lambda = lambda.res, Vorg = Vorg, V = V, 
      W = W, PI = PI, DELTA = DELTA, DELTA.bb = NULL, DELTA.ab = NULL, 
      DELTA.aa.b = NULL, GAMMA = GAMMA, test.name = "Johansen-Procedure")
}

H <- matrix(c(1,0,0,1,1), nrow = 5)

blrtest.jomoni(z = z, H = H, r = 1)
summary(blrtest.jomoni(z = z, H = H, r = 1))
