#' Random projection
#'
#' @param X numeric matrix of \code{n} by \code{d}
#' @param k numeric
#' @param type character
#'
#' @return numeric matrix
#' @export
random_projection <- function(X, k, type = "binary"){
  d <- ncol(X)
  stopifnot(d >= k, k > 0, k %% 1 == 0)

  if(type == "binary"){
    proj_mat <- .generate_binary_projection(d,k)
  } else if (type == "gaussian"){
    proj_mat <- .generate_gaussian_projection(d,k)
  } else stop()

  (X %*% proj_mat)/sqrt(k)
}

check_pairwise_distances <- function(X, Z, epsilon = 1e-4){
  stopifnot(nrow(X) == nrow(Z))

  dist_org <- as.matrix(stats::dist(X))^2
  dist_new <- as.matrix(stats::dist(Z))^2

  all((1-epsilon)*dist_org <= dist_new) &  all(dist_new <= (1+epsilon)*dist_org)
}


#########################

.generate_binary_projection <- function(d, k){
  val <- sample(c(-1,0,1), size = d*k, replace = T,
                prob = c(1/6,2/3,1/6))
  sqrt(3)*matrix(val, ncol = k)
}

.generate_gaussian_projection <- function(d, k){
  matrix(stats::rnorm(d*k), ncol = k)
}
