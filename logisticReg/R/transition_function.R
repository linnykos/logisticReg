# vec = sapply(seq(0, 20, length.out = 100), function(x){h_mle(0, x)})
# plot(vec, seq(0, 20, length.out = 100), xlim = c(0, 0.6), ylim = c(0, 20), asp = 0.6/20)

h_mle <- function(beta_0, gamma_0, n = 100000, tol = 1e-4, iter_max = 100){
  x <- stats::rnorm(n = n, mean = 0, sd = 1)
  val <- beta_0 + x*gamma_0
  prob_vec <- exp(val)/(1+exp(val))
  y <- sapply(1:n, function(i){sample(c(1,-1), 1, prob = c(prob_vec[i], 1-prob_vec[i]))})

  v <- y*x

  z <- stats::rnorm(n = n, mean = 0, sd = 1)

  t_0_prev <- Inf; t_1_prev <- Inf
  t_0 <- 0; t_1 <- 0
  iter <- 1
  while(iter < iter_max & abs(sum(c(t_0_prev, t_1_prev) - c(t_0, t_1))) > tol){
    t_0_prev <- t_0; t_1_prev <- t_1

    t_1 <- .alternate_minimization(y, v, z, t_0_init = t_0_prev, t_1_init = NA)
    t_0 <- .alternate_minimization(y, v, z, t_0_init = NA, t_1_init = t_1)

    c(t_0, t_1)

    iter <- iter+1
  }

  mean(pmax(t_0*y + t_1*v - z, 0)^2)
}

###############

.alternate_minimization <- function(y, v, z, t_0_init = NA, t_1_init = NA){
  stopifnot(!is.na(t_0_init) | !is.na(t_1_init))

  n <- length(y)

  if(is.na(t_1_init)){
    fn <- function(t_1){
      mean(pmax(t_0_init*y + t_1*v - z, 0)^2)
    }
  } else {
    fn <- function(t_0){
      mean(pmax(t_0*y + t_1_init*v - z, 0)^2)
    }
  }

  res <- stats::optim(0, fn = fn, lower = -500, upper = 500, method = "Brent")
  res$par
}
