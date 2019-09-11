rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 1

seq_val <- seq(0, 1, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)
y <- y_mat[116,]

# .distance_to_stability_set5 (dat, y, lambda, distr_class = "bernoulli")
###############
distr_class = "bernoulli"

stopifnot(length(y) == nrow(dat))
y <- as.numeric(y)
n <- nrow(dat); d <- ncol(dat)

powerset <- .powerset(1:d)[-1] # determine b_idx

all_dist_vec <- unlist(lapply(powerset, function(b_idx){
  print(paste0("b_idx: ", paste0(b_idx, collapse = "-")))
  len <- length(b_idx)
  powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
  k <- length(powerset_inner)

  dist_vec <- rep(NA, k^2)

  for(i in 1:k){
    for(j in 1:k){
      print(paste0("i: ", i))
      print(paste0("j: ", j))
      s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
      a_idx <- b_idx[powerset_inner[[j]]]

      kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
      plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)
      if(distr_class == "gaussian") {
        point <- .projection_euclidean(rep(0,n), plane)
        point <- .conjugate_grad_gaussian(point)

      } else if(distr_class == "bernoulli") {
        point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
        if(all(is.na(point)) | any(point <= 0) | any(point >= 1)) {
          dist_vec[(i-1)*k+j] <- NA; next()
        }
        point <- .conjugate_grad_bernoulli(point)

      } else {
        stop("distr_class not properly specified")
      }


      mab <- .construct_mab(dat, a_idx = a_idx, b_idx = b_idx)

      if(all(mab == 0) | length(mab) == 0){
        dist_vec[(i-1)*k+j] <- NA
      } else {
        nullspace <- .nullspace(mab)
        nullspace <- .plane(basis = nullspace)

        dist_vec[(i-1)*k+j] <- .distance_point_to_plane(point, nullspace)
      }
    }
  }

  dist_vec
}))
########################

b_idx = powerset[[1]]
i=2
j=1
  len <- length(b_idx)
  powerset_inner <- .powerset(1:len) # determine s_vec OR a_idx
  k <- length(powerset_inner)

s_vec <- rep(-1, len); s_vec[powerset_inner[[i]]] <- 1
a_idx <- b_idx[powerset_inner[[j]]]

kbs <- .construct_kbs(dat, b_idx, s_vec, lambda)
plane <- .plane(basis = kbs$basis, offset = y - kbs$offset)

# point <- .projection_bregman(rep(0.5,n), plane, distr_class = distr_class)
##########################

z = rep(0.5,n)
max_iter = 100
tol = 1e-3
if(nrow(plane$A) == ncol(plane$A)) {
  res <- as.numeric(plane$b)
  attributes(res) <- list("iteration" = NA, "tolerance" = NA)
  return(res)
}
res <- .setup_bregman_function(plane, distr_class)
f <- res$f; grad_f <- res$grad_f; prox <- res$prox
grad_f_z <- grad_f(z)

g <- function(x){f(x) - f(z) - t(grad_f_z) %*% (x-z)}
grad_g <- function(x){grad_f(x) - grad_f_z}
G_t <- function(x, eta){(x - prox(x - eta * (grad_f(x) - grad_f_z)))/eta}

iter <- 1
x_prev <- rep(Inf, length(z))
x_current <- res$x_current
if(all(is.na(x_current))) {
  res <- as.numeric(x_current)
  attributes(res) <- list("iteration" = NA, "tolerance" = NA)
  return(res)
}
