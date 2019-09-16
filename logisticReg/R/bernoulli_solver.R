bernoulli_solver <- function(dat, y, lambda, polytope = NA){
  stopifnot(all(dim(dat) == c(2,3)), length(y) == 2, all(y >= 0), all(y <= 1))

  if(all(is.na(polytope))){
    polytope <- .form_polytope(dat, lambda)

    # shift polytope
    for(i in 1:length(polytope)){
      tmp <- .plane(basis = polytope[[i]]$plane$basis,
                                    offset = polytope[[i]]$plane$offset + rep(0.5, 2))
      attr(tmp, "model") <- attr(polytope[[i]]$plane, "model")
      polytope[[i]]$plane <- tmp
      polytope[[i]]$intersection_1 <- polytope[[i]]$intersection_1 + rep(0.5, 2)
      polytope[[i]]$intersection_2 <- polytope[[i]]$intersection_2 + rep(0.5, 2)
    }
  }

  # determine if the point is inside the polytope
  bool_vec <- sapply(1:length(polytope), function(i){
    all(sign(polytope[[i]]$plane$A %*% rep(.5, 2) - polytope[[i]]$plane$b) ==
          sign(polytope[[i]]$plane$A %*% y - polytope[[i]]$plane$b))
  })
  if(all(bool_vec)){
    return(list(point = y, model = rep(NA, 3)))
  }

  # bregman project onto each of the faces of the polytope
  point_mat <- t(sapply(1:length(polytope), function(i){
    res <- .projection_bregman(y, polytope[[i]]$plane, distr_class = "bernoulli")
  }))

  # compute the bregman divergence at each point, include the corners of polytope
  for(i in 1:length(polytope)){
    point_mat <- rbind(point_mat, polytope[[i]]$intersection_1, polytope[[i]]$intersection_2)
  }
  val <- apply(point_mat, 1, function(x){
    .conjugate_bernoulli(x) - .conjugate_bernoulli(y) - (.conjugate_grad_bernoulli(y))%*%(x-y)
  })

  mat <- cbind(point_mat, val, 1:nrow(point_mat))
  colnames(mat) = c("x1", "x2", "val", "idx")

  # determine which index is appropriate
  while(TRUE){
    idx <- which.min(mat[,"val"])
    if(mat[idx, "idx"] > length(polytope)) break()

    polytope_idx <- mat[idx, "idx"]
    vec <- mat[idx, 1:2]
    if(sign(vec[1] - polytope[[polytope_idx]]$intersection_1[1]) != sign(vec[1] - polytope[[polytope_idx]]$intersection_2[1]) &
       sign(vec[2] - polytope[[polytope_idx]]$intersection_1[2]) != sign(vec[2] - polytope[[polytope_idx]]$intersection_2[2])) break()
    mat <- mat[-idx,]
  }

  # collect the model_vec
  if(idx <= length(polytope)){
    model_vec <- attr(polytope[[idx]]$plane, "model")
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
