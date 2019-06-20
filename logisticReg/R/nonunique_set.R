.construct_mab <- function(dat, a_idx, b_idx, tol = 1e-6){
  stopifnot(all(a_idx %in% b_idx), max(b_idx) <= ncol(dat))
  n <- nrow(dat); d <- ncol(dat)
  if(all(b_idx %in% a_idx)) return(matrix(0, ncol = 1, nrow = n)) # not sure about this

  minusb <- c(1:d)[-b_idx]
  d_minusb <- .diagonal_matrix(d, minusb)
  null_d_minusb <- .nullspace(d_minusb)

  d_bnota <- .diagonal_matrix(d, setdiff(b_idx, a_idx))
  # compute D_{B \backslash A} (X P_{null(D_{-B})})^+
  term2 <- d_bnota %*% MASS::ginv(dat %*% .projection_matrix(null_d_minusb))

  null_dat <- .nullspace(dat)
  idx <- sort(unique(unlist(lapply(1:ncol(null_d_minusb), function(x){
    which(abs(null_d_minusb[,x]) >= tol)
  }))))

  null_intersect <- null_dat
  if(length(idx) < nrow(null_dat)){
    null_intersect[-idx,] <- 0
  }
  null_intersect <- apply(null_intersect, 2, function(x){x/.l2norm(x)})
  projection_space <- .orthogonal_basis(d_bnota %*% null_intersect)
  stopifnot(nrow(projection_space) <= ncol(projection_space))

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
