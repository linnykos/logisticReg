#compute Lemma 6 of Tibs Taylor 2012: Degrees of freedom in lasso problems

.distance_to_stability_set2 <- function(dat, y, lambda){
  stopifnot(length(y) == nrow(dat))
  y <- as.numeric(y)
  n <- nrow(dat); d <- ncol(dat)

  powerset <- .powerset(1:d)[-1] # determine e_idx

  all_dist_vec <- unlist(lapply(powerset, function(e_idx){
    len <- length(e_idx)
    powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
    k <- length(powerset_inner)

    dist_vec <- rep(NA, k^2)

    for(i in 1:k){
      for(j in 1:k){
        s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
        a_idx <- powerset_inner[[j]] # here, a_idx represents the vector of positions locally within e_idx

        z_mat <- .compute_z2(dat, e_idx, a_idx)
        if(all(z_mat == 0)){
          dist_vec[(i-1)*k+j] <- NA
        } else {
          dist_vec[(i-1)*k+j] <- .l2norm(z_mat %*% (y - lambda*t(MASS::ginv(dat[,e_idx, drop = F]))%*%s_vec))
        }
      }
    }

    min(dist_vec, na.rm = T) # min over all possible sign vectors and a_idx
  }))

  min(all_dist_vec, na.rm = T) # min over all e_idx
}

.compute_z2 <- function(dat, e_idx, a_idx){
  stopifnot(length(a_idx) <= length(e_idx))
  if(length(a_idx) > 1) stopifnot(max(a_idx) <= length(e_idx))
  if(length(a_idx) == length(e_idx)) return(0)

  null <- .nullspace(dat[, e_idx, drop = F])

  if(length(a_idx) == 0){
    invert_idx <- 1:length(e_idx)
  } else {
    invert_idx <- (1:length(e_idx))[-a_idx]
  }

  proj_mat <- .projection_matrix(.orthogonal_basis(null[invert_idx, ,drop = F]))

  proj_mat %*% (MASS::ginv(dat[, e_idx, drop = F])[invert_idx, , drop = F])
}
