rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 1

seq_val <- seq(0, 1, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

dist_vec5 <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  .distance_to_stability_set5 (dat, as.numeric(y_mat[i,]), lambda, distr_class = "bernoulli")
})

bool_vec <- (dist_vec5 <= 1e-1)

plot(y_mat[,1], y_mat[,2], pch = 16, col = bool_vec+1, asp = T)
