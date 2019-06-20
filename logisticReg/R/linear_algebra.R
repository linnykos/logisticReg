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

#nullspace is the complement of the rowspace
#assumes no linearity in the columns
.nullspace <- function(mat){
  if(nrow(mat) == 0) return(diag(ncol(mat)))
  # if(ncol(mat) <= nrow(mat)) return(matrix(NA, nrow = ncol(mat), ncol = 0))

  .orthogonal_basis(t(mat))
}
