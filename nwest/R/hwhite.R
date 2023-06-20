hwhite <- function(y, x) {
  
  if (missing(y) | missing(x)) {
    stop("Missing arguments")
  }
  
  nobs <- nrow(x)
  nvar <- ncol(x)
  
  results <- list()
  results$meth <- 'hwhite'
  results$y <- y
  results$nobs <- nobs
  results$nvar <- nvar
  
  r <- qr.R(qr(x))
  xpxi <- solve(t(r) %*% r)
  
  results$beta <- xpxi %*% t(x) %*% y
  results$yhat <- x %*% results$beta
  results$resid <- y - results$yhat
  sigu <- t(results$resid) %*% results$resid
  results$sige <- sigu / (nobs - nvar)
  
  # perform White's correction
  xuux <- mcov(x, results$resid)
  xpxia <- xpxi %*% xuux %*% xpxi
  tmp <- sqrt(diag(xpxia))
  results$tstat <- results$beta / tmp
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

mcov <- function(x, u) {
  if (missing(x) | missing(u)) {
    stop("Wrong # of arguments to mcov")
  }
  
  nobs <- nrow(x)
  nvar <- ncol(x)
  
  xuux <- matrix(0, nvar, nvar)
  
  for (i in 1:nobs) {
    xp <- x[i, ]
    xpx <- xp %*% t(xp)
    upu <- u[i] * u[i]
    xuux <- xuux + upu * xpx
  }
  
  return(xuux)
}
