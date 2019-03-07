#' Determine if there exists a solution to logistic regression
#'
#' Returns \code{TRUE} if there exists a solution, \code{FALSE} if
#' not.
#'
#' @param X numeric matrix of \code{n} by \code{d}
#' @param y numeric vector of length \code{n}
#' @param tol positive numeric
#'
#' @return boolean
#' @export
#'
#' @examples
#' set.seed(10)
#' n <- 1000; d <- 10
#' X <- MASS::mvrnorm(n, rep(0,d), diag(d))
#' y <- generate_y_from_x(X)
#' existence(X, y) # should return TRUE
#'
#' set.seed(10)
#' n <- 1000; d <- 700
#' X <- MASS::mvrnorm(n, rep(0,d), diag(d))
#' beta <- rep(10/sqrt(d), d)
#' y <- generate_y_from_x(X, beta = beta)
#' existence(X, y) # should return FALSE
existence <- function(X, y, tol = 1e-3){
  stopifnot(nrow(X) == length(y))

  n <- nrow(X); d <- ncol(X)
  obj <- c(rep(0, d+1), 1)

  tmp <- sapply(1:n, function(i){y[i]*X[i,,drop = F]})
  if(ncol(X) > 1){
    tmp <- t(tmp)
  }

  A <- cbind(y, tmp, 1)
  res <- clplite::clp_solve(obj, A,
                            constr_lb = rep(tol, n),
                            constr_ub = rep(Inf, n),
                            var_lb = c(tol, rep(-Inf, d), 0),
                            var_ub = rep(Inf, d+2), max = F)

  # vec <- sapply(1:n, function(i){
  #   y[i]*(res$solution[1] + X[i,,drop=F]%*%res$solution[2:(d+1)])
  # })

  if(res$solution[d+2] < tol) FALSE else TRUE
}
