#' Generate Y from X
#'
#' According to logistic model
#'
#' @param X numeric matrix of \code{n} by \code{d}
#' @param beta_0 numeric
#' @param beta numeric vector of length \code{d}
#'
#' @return numeric vector for \code{c(-1,1)} of length \code{n}
#' @export
generate_y_from_x <- function(X, beta_0 = 0,
                              beta = rep(1/sqrt(ncol(X)), ncol(X))){
  n <- nrow(X)
  prob_vec <- sapply(1:n, function(i){
    val <- beta_0 + X[i,]%*%beta
    exp(val)/(1+exp(val))
  })

  sapply(1:n, function(i){sample(c(1,-1), 1, prob = c(prob_vec[i], 1-prob_vec[i]))})
}

