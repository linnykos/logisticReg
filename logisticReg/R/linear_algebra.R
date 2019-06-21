.l2norm <- function(x){sqrt(sum(x^2))}

.projection_matrix <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  if(ncol(mat) == 0) return(matrix(0, nrow(mat), nrow(mat)))

  mat %*% solve(t(mat) %*% mat) %*% t(mat)
}

.orthogonalize <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  d <- ncol(mat); n <- nrow(mat)

  mat[,1] <- mat[,1]/.l2norm(mat[,1])
  if(ncol(mat) == 1) return(mat)

  for(i in 2:d){
    proj_mat <- .projection_matrix(mat[,c(1:(i-1)), drop = F])
    tmp <- (diag(n) - proj_mat) %*% mat[,i]
    mat[,i] <- tmp/.l2norm(tmp)
  }

  mat
}

.orthogonal_basis <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  n <- nrow(mat); d <- ncol(mat)
  if(n == d) return(matrix(NA, nrow = n, ncol = 0))

  rand_mat <- matrix(rnorm(n*(n-d)), nrow = n, ncol = n-d)
  proj_mat <- .projection_matrix(mat)
  rand_mat <- (diag(n) - proj_mat) %*% rand_mat

  rand_mat <- apply(rand_mat, 2, function(x){x/.l2norm(x)})
  .orthogonalize(rand_mat)
}

.is_vector_in_basis <- function(vec, mat, tol = 1e-16){
  stopifnot(ncol(mat) <= nrow(mat), length(vec) == nrow(mat))
  d <- ncol(mat)

  complement_mat <- .orthogonal_basis(mat)
  full_mat <- cbind(mat, complement_mat)

  sol <- solve(full_mat, vec)
  sol <- sol/.l2norm(sol)

  sum(abs(sol[-c(1:d)])) <= tol
}

.diagonal_matrix <- function(n, include_idx = c(1:n)){
  diag(n)[include_idx,,drop = F]
}

#assumes no linearity in the columns
.nullspace <- function(mat){
  if(nrow(mat) == 0) return(diag(ncol(mat)))

  MASS::Null(t(mat))
}

# represents the plane as Ax = b, but you pass in the column vectors (basis matrix)
# and a point the plane passes through.
# the two representations are
#  {y s.t. there exists alpha s.t. y = offset + basis %*% alpha}
# or
#  {y s.t. A %*% y = b}
.plane <- function(basis, offset = rep(0, nrow(basis))){
  stopifnot(length(offset) == nrow(basis))

  if(ncol(basis) > 0) {
    basis <- apply(basis, 2, function(x){x/.l2norm(x)})
    A <- t(.orthogonalize(basis)); b <- A %*% offset
  } else {
    A <- diag(nrow(basis)); b <- offset
  }

  structure(list(basis = basis, offset = offset,
                 A = A, b = b), class = "plane")
}

.point_on_plane <- function(plane){
  if(nrow(plane$A) == 1){
    d <- length(plane$A)
    vec <- rep(0, d)
    idx <- which(plane$A != 0)
    stopifnot(length(idx) >= 1)
    idx <- idx[1]
    vec[-idx] <- 1
    vec[idx] <- as.numeric(plane$b - plane$A[-idx]%*%vec[-idx])/plane$A[idx]

    vec
  } else {
    k <- nrow(plane$A); n <- ncol(plane$A)
    mat <- matrix(0, ncol = 3*n, nrow = k+2*n)

    mat[1:k,1:n] <- plane$A
    mat[1:k,(n+1):(2*n)] <- -plane$A

    diag(mat[(k+1):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1):(k+n),1:n]) <- -1

    diag(mat[(k+1+n):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1+n):nrow(mat),(n+1):(2*n)]) <- -1

    vec <- c(plane$b, rep(0, 2*n))
    res <- lpSolve::lp(objective.in = c(rep(0, 2*n), rep(1, n)), const.mat = mat,
                       const.dir = c(rep("=", k), rep(">=", 2*n)),
                       const.rhs = vec)

    if(res$status == 2) {
      stop("LP to find point on plane failed")
    }
    res$solution[1:n] - res$solution[(n+1):(2*n)]
  }
}

.distance_point_to_plane <- function(point, plane){
  stopifnot(length(point) == length(plane$A))
  stopifnot(nrow(plane$A) == 1)

  x <- .point_on_plane(plane)
  .l2norm(plane$A%*%(point - x))/.l2norm(plane$A)
}
