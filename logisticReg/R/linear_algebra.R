.l2norm <- function(x){sqrt(sum(x^2))}

#' Construct projection matrix
#'
#' Construct a projection matrix that projects onto the column space of
#' \code{mat}. Requires \code{mat} to have more rows than columns (i.e., low-dimension).
#'
#' @param mat matrix
#'
#' @return matrix
#' @export
.projection_matrix <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  if(ncol(mat) == 0) return(matrix(0, nrow(mat), nrow(mat)))

  mat %*% solve(t(mat) %*% mat) %*% t(mat)
}

#' Orthonormalize a matrix
#'
#' Given a matrix \code{mat}, output a matrix that represents
#' the same column space as \code{mat} but each column now
#' is orthogonal with one another and has unit norm.
#'
#' Return a matrix of 0 column if \code{mat} is an empty matrix
#'
#' @param mat matrix
#'
#' @return matrix
#' @export
.orthogonalize <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  if(ncol(mat) == 0) return(mat)
  if(all(dim(mat) == 1)) return(matrix(1,1,1))
  d <- ncol(mat); n <- nrow(mat)

  mat[,1] <- mat[,1]/.l2norm(mat[,1])
  if(ncol(mat) == 1) return(mat)

  for(i in 2:d){
    proj_mat <- .projection_matrix(mat[,c(1:(i-1)), drop = F])
    tmp <- mat[,i] - proj_mat %*% mat[,i]
    mat[,i] <- tmp/.l2norm(tmp)
  }

  mat
}

#' Construct the orthogonal basis
#'
#' Given a matrix, construct a orthonormal basis (i.e., basis
#' where the columns are orthogonal to each other and has unit norm)
#' of the vector space orthogonal to the column space of \code{mat}.
#'
#' This function uses random number generation.
#'
#' @param mat matrix
#'
#' @return matrix
#' @export
.orthogonal_basis <- function(mat){
  stopifnot(ncol(mat) <= nrow(mat))
  n <- nrow(mat); d <- ncol(mat)
  if(n == d) return(matrix(NA, nrow = n, ncol = 0))

  rand_mat <- matrix(stats::rnorm(n*(n-d)), nrow = n, ncol = n-d)
  proj_mat <- .projection_matrix(mat)
  rand_mat <- (diag(n) - proj_mat) %*% rand_mat

  rand_mat <- apply(rand_mat, 2, function(x){x/.l2norm(x)})
  if(length(rand_mat) == 1) rand_mat <- matrix(1, 1, 1)
  .orthogonalize(rand_mat)
}

.diagonal_matrix <- function(n, include_idx = c(1:n)){
  diag(n)[include_idx,,drop = F]
}

#' Find the nullspace of a matrix.
#'
#' Construct an orthonormal basis of the null space of \code{mat},
#' where it contains all vectors \code{x} such that \code{mat %*% x} equal
#' 0, i.e., left null space.
#'
#' This function assumes no linearity in the columns.
#'
#' @param mat matrix
#'
#' @return matrix
#' @export
.nullspace <- function(mat){
  if(nrow(mat) == 0) return(diag(ncol(mat)))

  MASS::Null(t(mat))
}

#' Construct a plane
#'
#' We use two different representations of the same plane, both outputed in the
#' list. The inputs to this function are \code{basis} (the column vectors representing
#' the basis matrix) and point \code{offset} that the plane passes through.
#'
#' The two representations are: 1) all vectors \code{y} such that
#' there exists a vector \code{alpha} such that \code{y = offset + basis %*% alpha},
#' and 2) all vectors \code{y} such that \code{A%*%y = b}.
#'
#' If you want to input \code{A} and \code{b} into \code{.plane}, you must
#' explicitly set both \code{basis} and \code{offset} to be \code{NA}.
#'
#' @param basis matrix
#' @param offset vector
#' @param A matrix
#' @param b vector
#'
#' @return list
#' @export
.plane <- function(basis, offset = rep(0, nrow(basis)),
                   A = NA, b = NA){
  # construct plane using basis and offset
  if(all(is.na(A)) & all(is.na(b))){
    stopifnot(length(offset) == nrow(basis))

    if(ncol(basis) > 0) {
      basis <- apply(basis, 2, function(x){x/.l2norm(x)})
      A <- t(.orthogonal_basis(basis)); b <- A %*% offset
    } else {
      A <- diag(nrow(basis)); b <- offset
    }

  # construct plane using A and b
  } else {
    stopifnot(all(is.na(basis)), all(is.na(offset)))
    stopifnot(nrow(A) <= ncol(A), nrow(A) == length(b))

    basis <- .nullspace(A)
    point_on_plane <- .point_on_plane(A, b)
    proj_mat <- .projection_matrix(basis)
    offset  <- point_on_plane - proj_mat %*% point_on_plane
  }

  structure(list(basis = basis, offset = offset,
                 A = A, b = b), class = "plane")
}

#' Find a point on the plane
#'
#' Find a point \code{x} on the plane
#' represented by \code{A%*% x = b}
#' by minimizing the L1 norm of \code{x}
#'
#' @param A matrix
#' @param b vector
#'
#' @return vector
#' @export
#'
#' @examples
.point_on_plane <- function(A, b){
  if(nrow(A) == 1){
    d <- length(A)
    vec <- rep(0, d)
    idx <- which(A != 0)
    stopifnot(length(idx) >= 1)
    idx <- idx[1]
    vec[-idx] <- 1
    vec[idx] <- as.numeric(b - A[-idx]%*%vec[-idx])/A[idx]

    vec
  } else {
    k <- nrow(A); n <- ncol(A)
    mat <- matrix(0, ncol = 2*n, nrow = k)

    mat[1:k,1:n] <- A
    mat[1:k,(n+1):(2*n)] <- -A

    constr_vec <- b
    obj_vec <- rep(1, 2*n)

    res <- clplite::clp_solve(obj = obj_vec, A = mat,
                              constr_lb = constr_vec, constr_ub = constr_vec,
                              var_lb = rep(0, 2*n), var_ub = rep(Inf, 2*n),
                              max = FALSE)

    if(res$status != 0) {
      stop("LP to find point on plane failed")
    }
    res$solution[1:n] - res$solution[(n+1):(2*n)]
  }
}

#' Compute the distance from a point to plane
#'
#' @param point vector
#' @param plane object of class \code{plane}
#'
#' @return numeric
#' @export
.distance_point_to_plane <- function(point, plane){
  stopifnot(length(point) == ncol(plane$A))

  # x <- .point_on_plane(plane)
  # .l2norm(plane$A%*%(point - x))/.l2norm(plane$A)

  proj_mat <- .projection_matrix(plane$basis)
  closest_point <- proj_mat%*%(point - plane$offset) + plane$offset
  .l2norm(point - closest_point)
}

# from https://math.stackexchange.com/questions/25371/how-to-find-basis-for-intersection-of-two-vector-spaces-in-mathbbrn
.intersect_bases <- function(basis1, basis2){
  super_basis <- cbind(basis1, -basis2)
  null <- MASS::Null(t(super_basis))

  if(ncol(null) == 0){
    intersection <- matrix(NA, nrow = nrow(basis1), ncol = 0)
  } else {
    intersection <- apply(null, 2, function(x){
      basis1 %*% x[1:ncol(basis1)]
    })
  }

  .orthogonalize(intersection)
}

# specialized to only two 1d lines
.intersect_lines <- function(plane1, plane2){
  stopifnot(class(plane1) == "plane", class(plane2) == "plane",
    all(dim(plane1$A) == c(1,2)), all(dim(plane2$A) == c(1,2)),
    length(plane1$b) == 1, length(plane2$b) == 1)
  A_comb <- rbind(plane1$A, plane2$A)
  b_comb <- rbind(plane1$b, plane2$b)

  tryCatch({
    solve(A_comb, b_comb)
  }, error = function(e){
    rep(NA, 2)
  })
}
