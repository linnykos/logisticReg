rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

y <- c(0.85, 0.95)

###########################
polytope <- .form_polytope(dat, lambda)

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

point_mat <- t(sapply(1:length(polytope), function(i){
  res <- .projection_bregman(neg_grad_G, polytope[[i]]$plane, distr_class = "bernoulli",
                             offset = y, invert = T)
}))

# compute the bregman divergence at each point, include the corners of polytope
for(i in 1:length(polytope)){
  point_mat <- rbind(point_mat, polytope[[i]]$intersection_1, polytope[[i]]$intersection_2)
}

func1 <- .conjugate_bernoulli_constructor(offset = y, invert = T)
func2 <- .conjugate_grad_bernoulli_constructor(offset = y, invert = T)
val <- apply(point_mat, 1, function(x){
  if(all(is.na(x))) return(NA)
  func1(x) - func1(neg_grad_G) - func2(neg_grad_G)%*%(x-neg_grad_G)
})

mat <- cbind(point_mat, val, 1:nrow(point_mat))
colnames(mat) = c("x1", "x2", "val", "idx")

# determine which index is appropriate
while(TRUE){
  idx <- which.min(mat[,"val"])

  polytope_idx <- mat[idx, "idx"]
  if(mat[idx, "idx"] > length(polytope)) break()

  vec <- mat[idx, 1:2]
  if(sign(vec[1] - polytope[[polytope_idx]]$intersection_1[1]) != sign(vec[1] - polytope[[polytope_idx]]$intersection_2[1]) &
     sign(vec[2] - polytope[[polytope_idx]]$intersection_1[2]) != sign(vec[2] - polytope[[polytope_idx]]$intersection_2[2])) break()
  mat <- mat[-idx,]
}

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

list(point = mat[idx,1:2], model = model_vec)


# point_mat <- t(sapply(1:length(polytope), function(i){
#   res <- .projection_bregman(neg_grad_G, polytope[[i]]$plane, offset = y, distr_class = "bernoulli")
# }))
#
# f <- .conjugate_bernoulli_constructor(); grad_f <- .conjugate_grad_bernoulli_constructor()

# for(i in 1:6){
#   z = neg_grad_G
#   plane = polytope[[i]]$plane
#   offset = y
#   distr_class = "bernoulli"
#   max_iter = 100
#   tol = 1e-3
#   invert = T
#
#   if(nrow(plane$A) == ncol(plane$A)) {
#     res <- as.numeric(plane$b)
#     attributes(res) <- list("iteration" = NA, "tolerance" = NA)
#     return(res)
#   }
#   res <- .setup_bregman_function(plane, distr_class, offset = offset, invert = invert)
#   print(paste0(i, ": ", paste0(res$x_current, collapse = ", ")))
# }
#
#
#
# ################################
#
plot(NA, xlim = c(-1,1), ylim = c(-1,1), asp = T)
for(i in 1:length(polytope)){
  zz = rnorm(1000)*10
  val = sapply(1:length(zz), function(x){polytope[[i]]$plane$basis*zz[x]+polytope[[i]]$plane$offset})
  lines(val[1,], val[2,], col = i, lwd = 2)
}
points(neg_grad_G[1], neg_grad_G[2], col = "red", pch = 16)

model_vec
for(i in 1:6){
  print(paste0("Plane ", i, ": ", paste0(attr(polytope[[i]]$plane, "model"), collapse = ",")))
}
