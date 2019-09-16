rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

y <- c(0.85, 0.95)

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

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)

for(i in 1:length(polytope)){
  zz = rnorm(1000)*10
  val = sapply(1:length(zz), function(x){polytope[[i]]$plane$basis*zz[x]+polytope[[i]]$plane$offset})
  lines(val[1,], val[2,], col = i, lwd = 2)
}

#intersection is plane 1 (left: -1,0,0) and plane 6 (right: 0,0,1), intersection point is (0.5723, 0.7586)
#so we want to set b_idx = c(1,3) with s_vec = c(-1,1), and a_idx = 3
######################

distr_class = "bernoulli"

y <- as.numeric(y)
n <- nrow(dat); d <- ncol(dat)

powerset <- .powerset(1:d)[-1] # determine b_idx
b_idx <- powerset[[5]]

len <- length(b_idx)
powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
k <- length(powerset_inner)

i <- 3
j <- 3

s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
a_idx <- b_idx[powerset_inner[[j]]]

kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)

point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
point <- .conjugate_grad_bernoulli(point)
mab <- .construct_mab(dat, a_idx = a_idx, b_idx = b_idx)

nullspace <- .nullspace(mab)
nullspace <- .plane(basis = nullspace)

.distance_point_to_plane(point, nullspace)

###################

# we can hope to reverse this: what are all the points such that when you take the conjugate gradient,
## it lands in the nullspace

nullspace <- .nullspace(mab)
nullspace <- .plane(basis = nullspace)

tmp <- sort(rnorm(100000, sd = 100))
points_in <- sapply(tmp, function(x){
  zz <- nullspace$offset + x*nullspace$basis
  exp(zz)/(1+exp(zz))
})

plot(points_in[1,], points_in[2,], asp = T)
lines(rep(plane$offset[2], 2), c(0,1), col = "red")
lines(c(0,1), rep(plane$offset[1],2), col = "red")

########################3

rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

seq_val <- seq(0, 1, length.out = 100)
seq_val <- seq_val[-c(1, length(seq_val))]
y_mat <- expand.grid(seq_val, seq_val)

dist_vec5 <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  distr_class = "bernoulli"

  y <- as.numeric(y_mat[i,])
  n <- nrow(dat); d <- ncol(dat)

  powerset <- .powerset(1:d)[-1] # determine b_idx
  b_idx <- powerset[[5]]

  len <- length(b_idx)
  powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
  k <- length(powerset_inner)

  i <- 3
  j <- 3

  s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
  a_idx <- b_idx[powerset_inner[[j]]]

  kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
  plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)

  point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
  if(all(is.na(point)) | any(point <= 0) | any(point >= 1)) {
    return(NA)
  }
  point <- .conjugate_grad_bernoulli(point)
  mab <- .construct_mab(dat, a_idx = a_idx, b_idx = b_idx)

  nullspace <- .nullspace(mab)
  nullspace <- .plane(basis = nullspace)

  .distance_point_to_plane(point, nullspace)
})

bool_vec <- dist_vec5 <= 1e-1/3
plot(y_mat[,1], y_mat[,2], pch = 16, col = c("gray85", "red")[bool_vec+1], asp = T)

