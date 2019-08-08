#compute Lemma 5 of Tibs Taylor 2012: Degrees of freedom in lasso problems

.distance_to_stability_set1 <- function(dat, y, lambda){
  stopifnot(length(y) == nrow(dat))
  y <- as.numeric(y)
  n <- nrow(dat); d <- ncol(dat)

  powerset <- .powerset(1:d)[-1] # determine e_idx

  all_dist_vec <- unlist(lapply(powerset, function(e_idx){
    len <- length(e_idx)
    powerset_inner <- .powerset(1:len) # determine s_vec
    k <- length(powerset_inner)

    dist_vec <- rep(NA, k)

    for(i in 1:k){
      s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
      dist_vec[i] <- .distance_to_stability_set1_instance(dat, y, lambda, e_idx, s_vec)
    }

    min(dist_vec) # min over all possible sign vectors
  }))

  min(all_dist_vec) # min over all e_idx
}

.distance_to_stability_set1_instance <- function(dat, y, lambda, e_idx, s_vec){
  ginv_mat <- MASS::ginv(dat[,e_idx,drop = F])

  min(abs(sapply(1:len, function(x){
    ginv_mat[x,] %*% (y - lambda*t(ginv_mat)%*%s_vec)
  })))
}
