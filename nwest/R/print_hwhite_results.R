print_hwhite_results <- function(results, vnames) {
  cat('\n')
  cat('White Heteroscedastic Consistent Estimates\n')
  
  if (!is.null(vnames)) {
    cat('Dependent Variable: ', vnames[1], '\n')
  }
  
  cat('R-squared: ', round(results$rsqr, 4), '\n')
  cat('R-bar-squared: ', round(results$rbar, 4), '\n')
  cat('Durbin-Watson Statistic: ', round(results$dw, 4), '\n')
  cat('Number of observations, variables: ', results$nobs, ', ', results$nvar, '\n')
  
  cat('\n')
  cat('Coefficients, T-statistics, T-probabilities:\n')
  coef_names <- vnames[2:length(vnames)]
  
  max_name_length <- max(nchar(coef_names))
  
  for (i in seq_along(coef_names)) {
    coef_str <- sprintf('%*.4f', max_name_length, results$beta[i])
    tstat_str <- sprintf('%*.4f', max_name_length, results$tstat[i])
    tprob_str <- sprintf('%*.4f', max_name_length, pt(abs(results$tstat[i]), results$nobs - results$nvar, lower.tail = FALSE) * 2)
    cat(paste0(sprintf('%-*s', max_name_length, coef_names[i]), ': Coefficient: ', coef_str, ', T-statistic: ', tstat_str, ', T-probability: ', tprob_str, '\n'))
  }
  
  cat('***************************************************************\n')
}
