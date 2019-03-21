#' Computing the margin
#'
#' @param X numeric matrix of \code{n} by \code{d}
#' @param y numeric vector of length \code{n}
#' @param check boolean
#' @param tol small positive numeric
#'
#' @return numeric
#' @export
#'
#' @examples
#' set.seed(10)
#' n <- 1000; d <- 700
#' X <- MASS::mvrnorm(n, rep(0,d), diag(d))
#' beta <- rep(10/sqrt(d), d)
#' y <- generate_y_from_x(X, beta = beta)
#' existence(X, y) # should return FALSE, meaning there is a margin
#' margin(X, y)
margin <- function(X, y, check = TRUE, tol = 1e-2){
  res <- e1071::svm(x = X, y = y, type = "C-classification", cost = 1000000000,
                    scale = FALSE, kernel = "linear", shrinking = FALSE)
  if(check) stopifnot(all(y %in% c(-1,1)))

  beta = rep(0, 2)
  for(i in 1:length(res$index)){
    beta <- beta + res$coefs[i]*X[res$index[i],]
  }

  if(check){
    stopifnot(abs(sum(res$coefs)) < tol,
              length(res$coefs) == length(res$index))

    #checking the distance from margin is constant
    # it's unclear to me if the labels are sometimes flipped...
    beta_0_vec1 <- sapply(1:length(res$index), function(i){
      X[res$index[i],]%*%beta+y[res$index[i]]
    })
    beta_0_vec2 <- sapply(1:length(res$index), function(i){
      X[res$index[i],]%*%beta-y[res$index[i]]
    })

    stopifnot((max(beta_0_vec1) - min(beta_0_vec1) < tol) |
                (max(beta_0_vec2) - min(beta_0_vec2) < tol))
  }

  2/.l2norm(beta)
}

###############

.l2norm <- function(x){sqrt(sum(x^2))}
