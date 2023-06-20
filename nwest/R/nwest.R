nwest <- function(y, x, nlag) {
  
  if (missing(y) | missing(x) | missing(nlag)) {
    stop("Missing arguments")
  }
  
  nobs <- nrow(x)
  nvar <- ncol(x)
  
  results <- list()
  results$meth <- 'nwest'
  results$y <- y
  results$nobs <- nobs
  results$nvar <- nvar
  
  xpxi <- solve(t(x) %*% x)
  results$beta <- xpxi %*% t(x) %*% y
  results$yhat <- x %*% results$beta
  results$resid <- y - results$yhat
  sigu <- t(results$resid) %*% results$resid
  results$sige <- sigu/(nobs - nvar)
  
  # perform Newey-West correction
  emat <- matrix(nrow = nvar, ncol = nobs)
  for (i in 1:nvar){
    emat[i, ] <- results$resid
  }
  
  hhat <- emat * t(x)
  G <- matrix(0, nrow = nvar, ncol = nvar)
  w <- rep(0, (2 * nlag) + 1)
  a <- 0

  while (a != nlag + 1) {
    # initialize ga as a zero matrix
    ga <- matrix(0, nrow = nvar, ncol = nvar)
    w[nlag + 1 + a] <- (nlag + 1 - a) / (nlag + 1)
    za <- hhat[ , (a + 1):nobs] %*% t(hhat[ , 1:(nobs - a)])
    if (a == 0){
      ga <- ga + za
    } else {
      ga <- ga + za + t(za)
    }
    G <- G + w[nlag + 1 + a] * ga
    a <- a + 1
  }
  
  V <- xpxi %*% G %*% xpxi
  nwerr <- sqrt(diag(V))
  
  results$tstat <- results$beta / nwerr # Newey-West t-statistics
  ym <- y - mean(y)
  rsqr1 <- sigu
  rsqr2 <- t(ym) %*% ym
  results$rsqr <- 1.0 - rsqr1 / rsqr2 # r-squared
  rsqr1 <- rsqr1 / (nobs - nvar)
  rsqr2 <- rsqr2 / (nobs - 1.0)
  results$rbar <- 1 - (rsqr1 / rsqr2) # rbar-squared
  ediff <- results$resid[2:nobs] - results$resid[1:(nobs - 1)]
  sigu_vec <- rep(sigu, length(ediff))
  results$dw <- sum((ediff^2) / sigu_vec) # durbin-watson
  
  return(results)
}
