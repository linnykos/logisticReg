rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.5

###############3

stopifnot(all(dim(dat) == c(2,3)))

n <- nrow(dat); d <- ncol(dat)

# populate the polytope
point_mat <- .populate_polytope(dat, lambda)
plot(point_mat[,1], point_mat[,2], asp = T)

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
    counter <- counter + 1
  }
}

for(i in 1:length(kbs_list)){
  zz = rnorm(1000)*10
  val = sapply(1:length(zz), function(x){kbs_list[[i]]$basis*zz[x]+kbs_list[[i]]$offset})
  lines(val[1,], val[2,], col = "green", lwd = 2)
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

for(i in 1:6){
  points(lis[[i]]$intersection_1[1], lis[[i]]$intersection_1[2], col = "red", pch = 16)
  points(lis[[i]]$intersection_2[1], lis[[i]]$intersection_2[2], col = "red", pch = 16)
}
