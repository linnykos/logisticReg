.projection_euclidean <- function(point, plane){
  stopifnot(class(plane) == "plane")
  if(ncol(plane$basis) == 0) return(plane$offset)

  point2 <- point - plane$offset
  proj_mat <- .projection_matrix(plane$basis)
  proj_point <- proj_mat %*% point2

  proj_point + plane$offset
}
