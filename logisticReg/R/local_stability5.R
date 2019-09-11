# computing Remark 10 in Ali, Tibs 2018: The generalized lasso problem and uniqueness

.distance_to_stability_set5 <- function(dat, y, lambda, distr_class = "gaussian"){
  stopifnot(length(y) == nrow(dat))
  y <- as.numeric(y)
  n <- nrow(dat); d <- ncol(dat)

  powerset <- .powerset(1:d)[-1] # determine b_idx

  all_dist_vec <- unlist(lapply(powerset, function(b_idx){
    # print(paste0("b_idx: ", paste0(b_idx, collapse = "-")))
    len <- length(b_idx)
    powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
    k <- length(powerset_inner)

    dist_vec <- rep(NA, k^2)

    for(i in 1:k){
      for(j in 1:k){
        # print(paste0("i: ", i))
        # print(paste0("j: ", j))
        s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
        a_idx <- b_idx[powerset_inner[[j]]]

        kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
        plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)
        if(distr_class == "gaussian") {
          point <- .projection_euclidean(rep(0,n), plane)
          point <- .conjugate_grad_gaussian(point)

        } else if(distr_class == "bernoulli") {
          point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
          if(all(is.na(point)) | any(point <= 0) | any(point >= 1)) {
            dist_vec[(i-1)*k+j] <- NA; next()
          }
          point <- .conjugate_grad_bernoulli(point)

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
}

.construct_mab <- function(dat, a_idx, b_idx, tol = 1e-6){
  stopifnot(all(a_idx %in% b_idx), max(b_idx) <= ncol(dat))
  n <- nrow(dat); d <- ncol(dat)
  if(all(b_idx %in% a_idx)) return(matrix(0, ncol = 1, nrow = n))

  minusb <- c(1:d)[-b_idx]
  d_minusb <- .diagonal_matrix(d, minusb)
  null_d_minusb <- .nullspace(d_minusb)

  d_bnota <- .diagonal_matrix(d, setdiff(b_idx, a_idx))
  # compute D_{B \backslash A} (X P_{null(D_{-B})})^+
  term2 <- d_bnota %*% MASS::ginv(dat %*% .projection_matrix(null_d_minusb))

  null_dat <- .nullspace(dat)
  null_intersect <- .intersect_bases(null_d_minusb, null_dat)

  projection_space <- .orthogonal_basis(d_bnota %*% null_intersect)

  # compute P_{[D_{B \backslash A}(null(X) \cap null(D_{-B}))]^c}
  term1 <- .projection_matrix(projection_space)

  term1 %*% term2
}

.construct_kbs <- function(dat, b_idx, s_vec, lambda){
  stopifnot(all(s_vec %in% c(-1,0,1)), length(s_vec) == length(b_idx))
  d <- ncol(dat)

  minusb <- c(1:d)[-b_idx]
  d_minusb <- .diagonal_matrix(d, minusb)
  d_b <- .diagonal_matrix(d, b_idx)

  proj_mat <- .projection_matrix(.nullspace(d_minusb))
  point <- lambda * MASS::ginv(proj_mat %*% t(dat)) %*% t(d_b) %*% s_vec

  plane <- .nullspace(proj_mat %*% t(dat))

  .plane(basis = plane, offset = point)
}

################

.powerset <- function(vec){
  if(length(vec) == 0) return(numeric(0))

  rje::powerSet(vec)
}
