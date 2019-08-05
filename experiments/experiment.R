rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 3

seq_val <- seq(-5, 5, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

############

combn_mat <- utils::combn(c(1:3), 2)
i <- 2
b_idx <- combn_mat[,i]
sign_combn <- .powerset(c(1,2))
j <- 2
s_vec <- rep(-1, 2)
s_vec[sign_combn[[j]]] <- 1

res <- .construct_kbs(dat, b_idx, s_vec, lambda)
res$offset

y <- res$offset
.distance_to_stability_set1 (dat, y, lambda)
