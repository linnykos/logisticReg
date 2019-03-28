#' Check gaussianity via marginals
#'
#' @param X numeric matrix of \code{n} by \code{d}
#' @param prob_vec numeric vector
#'
#' @return numeric vector
#' @export
gaussian_check <- function(X, prob_vec = c(0.5, 1)){
  d <- nrow(X)

  diff_vec <- apply(X, 2, function(x){
    stats::ks.test(x, stats::pnorm, alternative = "two.sided")$statistic
  })

  quantile(diff_vec, probs = prob_vec)
}
