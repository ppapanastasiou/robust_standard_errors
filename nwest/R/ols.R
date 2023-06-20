ols <- function(y, x) {
  if (length(dim(x)) != 2) {
    stop('x must be a matrix in ols')
  }
  nobs <- nrow(x)
  nvar <- ncol(x)
  
  if (length(y) != nobs) {
    stop('x and y must have the same number of observations in ols')
  }
  
  results <- list()
  results$meth <- 'ols'
  results$y <- y
  results$nobs <- nobs
  results$nvar <- nvar
  
  xpxi <- solve(crossprod(x))
  
  results$beta <- xpxi %*% t(x) %*% y
  results$yhat <- x %*% results$beta
  results$resid <- y - results$yhat
  sigu <- sum(results$resid^2)
  results$sige <- sigu / (nobs - nvar)
  tmp <- results$sige * diag(xpxi)
  sigb <- sqrt(tmp)
  results$bstd <- sigb
  results$tstat <- results$beta / sqrt(tmp)
  ym <- y - mean(y)
  rsqr1 <- sigu
  rsqr2 <- sum(ym^2)
  results$rsqr <- 1 - rsqr1 / rsqr2
  rsqr1 <- rsqr1 / (nobs - nvar)
  rsqr2 <- rsqr2 / (nobs - 1)
  
  if (rsqr2 != 0) {
    results$rbar <- 1 - rsqr1 / rsqr2
  } else {
    results$rbar <- results$rsqr
  }
  
  ediff <- results$resid[-1] - results$resid[-nobs]
  results$dw <- sum(ediff^2) / sigu
  
  return(results)
}
