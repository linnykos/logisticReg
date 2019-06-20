rm(list=ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 1.7
b_idx <- c(1,3)
s_vec <- c(-1,1)

kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
point <- .projection_euclidean(rep(0,2), .plane(basis = kbs$basis, offset = c(-10,10) - kbs$offset))

nullspace <- .nullspace(.construct_mab(dat, a_idx = 1, b_idx))

.is_vector_in_basis(point, nullspace)
