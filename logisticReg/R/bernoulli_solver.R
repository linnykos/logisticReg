bernoulli_solver <- function(dat, y, lambda, polytope = NA){
  stopifnot(all(dim(dat) == c(2,3)), length(y) == 2, all(y >= 0), all(y <= 1))

  if(all(is.na(polytope))){
    polytope <- .form_polytope(dat, lambda)
  }

  # determine if the point is inside the polytope
  bool_vec <- sapply(1:length(polytope), function(i){
    all(sign(polytope[[i]]$plane$A %*% rep(0, 2) - polytope[[i]]$plane$b) ==
          sign(polytope[[i]]$plane$A %*% (y - rep(0.5,2)) - polytope[[i]]$plane$b))
  })
  if(all(bool_vec)){
    return(list(point = y, model = rep(NA, 3)))
  }

  #compute -nalba G(0)
  neg_grad_G <- y - rep(0.5, 2)

  # bregman project onto each of the faces of the polytope
  point_mat <- t(sapply(1:length(polytope), function(i){
    .projection_bregman(neg_grad_G, polytope[[i]]$plane, distr_class = "bernoulli",
                        offset = y, invert = T,
                        intersection_1 = polytope[[i]]$intersection_1,
                        intersection_2 = polytope[[i]]$intersection_2)
  }))

  # compute the bregman divergence at each point, include the corners of polytope
  # put the intersection points above first, so the "which.min" function below hits them first
  for(i in 1:length(polytope)){
    point_mat <- rbind(polytope[[i]]$intersection_1, polytope[[i]]$intersection_2, point_mat)
  }
  func1 <- .conjugate_bernoulli_constructor(offset = y, invert = T)
  func2 <- .conjugate_grad_bernoulli_constructor(offset = y, invert = T)
  val <- apply(point_mat, 1, function(x){
    if(all(is.na(x))) return(NA)
    func1(x) - func1(neg_grad_G) - func2(neg_grad_G)%*%(x-neg_grad_G)
  })

  mat <- cbind(point_mat, val, c((length(polytope)+1):nrow(point_mat), 1:length(polytope)))
  colnames(mat) = c("x1", "x2", "val", "idx")

  # determine which index is appropriate
  idx <- which.min(mat[,"val"])
  polytope_idx <- mat[idx, "idx"]

  # collect the model_vec
  if(polytope_idx <= length(polytope)){
    model_vec <- attr(polytope[[polytope_idx]]$plane, "model")
  } else {
    plane_idx <- which(sapply(1:length(polytope), function(i){
      sum(abs(mat[idx,1:2] - polytope[[i]]$intersection_1)) == 0 | sum(abs(mat[idx,1:2] - polytope[[i]]$intersection_2)) == 0
    }))
    stopifnot(length(plane_idx) == 2)

    model_vec <- rep(NA, 3)
    for(i in 1:length(plane_idx)){
      tmp <- attr(polytope[[plane_idx[i]]]$plane, "model")
      model_vec[which(!is.na(tmp))] <- tmp[which(!is.na(tmp))]
    }
  }

  return(list(point = mat[idx,1:2], model = model_vec))
}
