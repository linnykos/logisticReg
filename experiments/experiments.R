rm(list=ls())
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 1.7
b_idx <- c(1,3)
s_vec <- c(-1,1)

kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
