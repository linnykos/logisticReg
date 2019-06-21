rm(list=ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 1.7
y <- c(5,5)

# .distance_point_to_set(dat, y, lambda)

############
n <- nrow(dat); d <- ncol(dat)

powerset <- .powerset(1:d)[-1]  # determine b_idx

b_idx <- powerset[[3]]
len <- length(b_idx)
powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
k <- length(powerset_inner)

dist_vec <- rep(NA, k^2)

for(i in 1:k){
  for(j in 1:k){
    s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
    a_idx <- b_idx[powerset_inner[[j]]]

    kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
    point <- .projection_euclidean(rep(0,n), .plane(basis = kbs$basis, offset = y - kbs$offset))

    mab <- .construct_mab(dat, a_idx = a_idx, b_idx = b_idx)

    if(all(mab == 0) | length(mab) == 0){
      dist_vec[(i-1)*k+j] <- NA
    } else {
      nullspace <- .nullspace(mab)
      nullspace <- .plane(basis = nullspace)

      dist_vec[(i-1)*k+j] <- .distance_point_to_plane(point, nullspace)
    }
  }
}
