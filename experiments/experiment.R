rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

seq_val <- seq(0, 1, length.out = 10)
y_mat <- expand.grid(seq_val, seq_val)

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

y <- c(0,0)

# bernoulli_solver (dat, y, lambda, polytope = polytope)
########################

# determine if the point is inside the polytope
bool_vec <- sapply(1:length(polytope), function(i){
  all(sign(polytope[[i]]$plane$A %*% rep(.5, 2) - polytope[[i]]$plane$b) ==
        sign(polytope[[i]]$plane$A %*% y - polytope[[i]]$plane$b))
})
if(all(bool_vec)){
  return(list(point = y, model = rep(NA, 3)))
}

# point_mat <- t(sapply(1:length(polytope), function(i){
#   print(i)
#   res <- .projection_bregman(y, polytope[[i]]$plane, distr_class = "bernoulli")
# }))

# .projection_bregman(y, polytope[[1]]$plane, distr_class = "bernoulli")
######################

z <- y
plane <- polytope[[1]]$plane
distr_class = "bernoulli"
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

# while(iter < max_iter & (is.na(tol) || .l2norm(x_prev - x_current) > tol)){
#   x_prev <- x_current
#
#   eta <- .backtrack_line_search(x_current, g, grad_g, G_t)
#   x_current <- prox(x_prev - eta * grad_g(x_prev))
#
#   iter <- iter + 1
#   # print(x_current)
# }

x_prev <- x_current

# eta <- .backtrack_line_search(x_current, g, grad_g, G_t)
######################
x <- x_current
beta = 0.5
eta_init = 1
tol = 1e-3

eta <- eta_init
gx <- g(x)
grad_gx <- grad_g(x)
Gtx <- G_t(x, eta)
# counter <- 1

# while(TRUE){
#   val1 <- g(x - eta*Gtx)
#   val2 <- gx - eta * t(grad_gx)%*%Gtx + eta*.l2norm(Gtx)^2/2
#   if(val2 > val1) break()
#   if(abs(eta) <= tol) {eta <- 0; break()}
#   # print(paste0(counter, ": ", eta, " // ", val1, " vs. ", val2))
#   eta <- eta*beta
#
#   counter <- counter + 1
# }

