# the following procedures are NOT general at all, and require a lot of hand-tuning
# this stems from my inability to represent the polytope in a convenient form that can be
#  automated (in terms of calculation)

.populate_polytope <- function(dat, lambda, grid_size = 100, lim = c(-.5,.5)){
  stopifnot(all(dim(dat) == c(2,3)))

  point_mat <- expand.grid(seq(lim[1], lim[2], length.out = grid_size), seq(lim[1], lim[2], length.out = grid_size))
  bool_vec <- apply(point_mat, 1, function(x){
    all(abs(t(dat) %*% x) <= lambda)
  })

  as.matrix(point_mat[bool_vec,])
}

# we use a really janky strategy: populate the polytope with a lot of points and then
#  form all possible k_bs lines. if a k_bs line has all the points in the polytope on one side,
#  then we say it's part of the polytopes
# to see if two lines are neighbors, see if their intersection yields a point near the polytope
.form_polytope <- function(dat, lambda){
  stopifnot(all(dim(dat) == c(2,3)))

  n <- nrow(dat); d <- ncol(dat)

  # populate the polytope
  point_mat <- .populate_polytope(dat, lambda)

  # form all possible kbs
  powerset <- list(1,2,3) # determine b_idx
  counter <- 1
  kbs_list <- vector("list", 1)

  for(b_idx in powerset){
    len <- length(b_idx)
    powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
    k <- length(powerset_inner)

    for(i in 1:k){
      s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
      kbs_list[[counter]] <- .construct_kbs(dat, b_idx, s_vec, lambda)
      model_vec <- rep(NA, 3); model_vec[b_idx] <- s_vec
      attr(kbs_list[[counter]], "model") <- model_vec
      counter <- counter + 1
    }
  }

  # determine which ones are on the polytope (they all should be)
  bool_vec <- sapply(1:length(kbs_list), function(i){
    vec <- point_mat %*% t(kbs_list[[i]]$A) - as.numeric(kbs_list[[i]]$b)
    all(sign(vec) == sign(vec[1]))
  })
  kbs_list <- kbs_list[which(bool_vec)]

  # determine the intersection of planes
  combn_list <- utils::combn(length(kbs_list),2)
  intersection_points <- apply(combn_list, 2, function(x){
    .intersect_lines(kbs_list[[x[1]]], kbs_list[[x[2]]])
  })

  # determine which interesctions are the corners of the polytope
  min_dist <- apply(intersection_points, 2, function(x){
    min(apply(point_mat, 1, function(y){.l2norm(y-x)}))
  })

  # determine which points are the intersections
  lis <- vector("list", length(kbs_list))
  for(i in 1:length(lis)){
    idx <- which(apply(combn_list, 2, function(x){i %in% x}))
    sorted_order <- order(min_dist[idx], decreasing = F, na.last = T)

    lis[[i]] <- list("plane" = kbs_list[[i]], intersection_1 = intersection_points[,idx[sorted_order[1]]],
                     intersection_2 = intersection_points[,idx[sorted_order[2]]])
  }

  structure(lis, class = "polytope")
}

