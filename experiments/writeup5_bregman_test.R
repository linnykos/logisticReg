rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 3

seq_val <- seq(0, 1, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

distr_class = "bernoulli"

##############

y <- rep(0.5, 2)
y <- as.numeric(y)
stopifnot(length(y) == nrow(dat))
y <- as.numeric(y)
n <- nrow(dat); d <- ncol(dat)

powerset <- .powerset(1:d)[-1] # determine b_idx

all_dist_vec <- unlist(lapply(powerset, function(b_idx){
  len <- length(b_idx)
  powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
  k <- length(powerset_inner)

  dist_vec <- rep(NA, k^2)

  for(i in 1:k){
    for(j in 1:k){
      s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
      a_idx <- b_idx[powerset_inner[[j]]]

      kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
      plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)
      if(distr_class == "gaussian") {
        point <- .projection_euclidean(rep(0,n), plane)
      } else if(distr_class == "bernoulli") {
        point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
        if(all(is.na(point))) {
          dist_vec[(i-1)*k+j] <- NA; next()
        }
      } else {
        stop("distr_class not properly specified")
      }

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

  dist_vec
}))

min(all_dist_vec, na.rm = T)
