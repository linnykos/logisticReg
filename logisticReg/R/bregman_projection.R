.plane <- function(basis, offset = rep(0, nrow(basis))){
  stopifnot(length(offset) == nrow(basis))

  basis <- apply(basis, 2, function(x){x/.l2norm(x)})

  structure(list(basis = basis, offset = offset), class = "plane")
}

.projection_euclidean <- function(point, plane){
  stopifnot(class(plane) == "plane")

  point2 <- point - plane$offset
  proj_mat <- .projection_matrix(plane$basis)
  proj_point <- proj_mat %*% point2

  proj_point + plane$offset
}
