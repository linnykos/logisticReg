set.seed(10)
cov_mat <- stats::toeplitz(c(5:1))
n <- 1000
X <- MASS::mvrnorm(n, rep(0, 5), cov_mat)
beta_0 <- 0.5
beta <- rep(1,5)
y <- apply(X, 1, function(x){
  p <- as.numeric(1/(1+exp(-beta_0 + x %*% beta)))
  sample(c(1,-1), 1, prob = c(p,1-p))
})

