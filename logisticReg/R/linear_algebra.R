.projection_matrix <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  mat %*% solve(t(mat) %*% mat) %*% t(mat)
}

.orthogonal_basis <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  n <- nrow(mat); d <- ncol(mat)
  if(n == d) return(matrix(NA, nrow = n, ncol = 0))

  rand_mat <- matrix(rnorm(n*(n-d)), nrow = n, ncol = n-d)
  proj_mat <- .projection_matrix(mat)
  rand_mat <- (diag(n) - proj_mat) %*% rand_mat

  apply(rand_mat, 2, function(x){x/.l2norm(x)})
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
  diag(n)[include_idx,]
}


