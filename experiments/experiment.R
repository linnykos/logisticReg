rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

y <- c(0.5, 0.76)

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

i=1
z = neg_grad_G
plane = polytope[[i]]$plane
distr_class = "bernoulli"
offset = y
invert = T
max_iter = 100
tol = 1e-3
intersection_1 <- polytope[[i]]$intersection_1
intersection_2 <- polytope[[i]]$intersection_2

if(nrow(plane$A) == ncol(plane$A)) {
  res <- as.numeric(plane$b)
  attributes(res) <- list("iteration" = NA, "tolerance" = NA)
  return(res)
}
res <- .setup_bregman_function(plane, distr_class, offset = offset, invert = invert)
x_current <- res$x_current

# f <- res$f; grad_f <- res$grad_f; prox <- res$prox
# grad_f_z <- grad_f(z)
#
# g <- function(x){f(x) - f(z) - t(grad_f_z) %*% (x-z)}
# grad_g <- function(x){grad_f(x) - grad_f_z}
# G_t <- function(x, eta){(x - prox(x - eta * (grad_f(x) - grad_f_z)))/eta}
#
# iter <- 1
# x_prev <- rep(Inf, length(z))
#
# x_prev <- x_current
#
# eta <- .backtrack_line_search(x_current, g, grad_g, G_t)
# x_current <- prox(x_prev - eta * grad_g(x_prev))
#
# iter <- iter + 1

#####################################

# beta = 0.5
# eta_init = 1
# tol = 1e-3
# x = x_current
#
# eta <- eta_init
# gx <- g(x)
# grad_gx <- grad_g(x)
# counter <- 1
#
# while(TRUE){
#   Gtx <- G_t(x, eta)
#   val1 <- g(x - eta*Gtx)
#   val2 <- gx - eta * t(grad_gx)%*%Gtx + eta*.l2norm(Gtx)^2/2
#   if(val2 > val1) break()
#   if(abs(eta) <= tol) {eta <- 0; break()}
#   print(paste0(counter, ": ", eta, " // ", val1, " vs. ", val2))
#   eta <- eta*beta
#
#   counter <- counter + 1
# }

######################

# try a simple experiment: create a grid along all feasible points and compute the bregman divergence
# can only handle the case where the plane is 1d in a 2d space
stopifnot(all(dim(plane$basis) == c(2,1)))
f <- .conjugate_bernoulli_constructor(offset = offset, invert = invert)
grad_f <- .conjugate_grad_bernoulli_constructor(offset = offset, invert = invert) #NOTE: y-1/2 (equal to offset-1/2 here) should be a valid input
domain_mat <- attr(grad_f, "domain")

plane_list <- list(.plane(basis = matrix(c(0,1), 2, 1), offset = c(domain_mat[1,1], domain_mat[2,1])),
                   .plane(basis = matrix(c(1,0), 2, 1), offset = c(domain_mat[1,1], domain_mat[2,1])),
                   .plane(basis = matrix(c(0,1), 2, 1), offset = c(domain_mat[1,2], domain_mat[2,2])),
                   .plane(basis = matrix(c(1,0), 2, 1), offset = c(domain_mat[1,2], domain_mat[2,2])))
point_mat <- sapply(plane_list, function(p){
  .intersect_two_lines(plane$A, p$A, plane$b, p$b)
})

idx <- which(apply(point_mat, 2, function(x){all(x >= domain_mat[,1]) & all(x <= domain_mat[,2])}))
point_mat <- point_mat[,idx]

# first check: make sure the point_mat is larger than the intersections
polytope[[1]]
# make sure the intersection points are along this line

zz <- polytope[[1]]$intersection_2
val1 <- (zz[1]-point_mat[1,1])/(point_mat[1,2] - point_mat[1,1])
recon1 <- point_mat[2,1] + val1*(point_mat[2,2] - point_mat[2,1])
sum(abs(recon1 - zz[2]))
#good

# make a grid of points, and compute the divergence at each point
x_seq <- seq(0, 1, length.out = 1000)
zz <- sapply(x_seq, function(x){point_mat[,1] + x * apply(point_mat, 1, diff)})

func1 <- .conjugate_bernoulli_constructor(offset = y, invert = T)
func2 <- .conjugate_grad_bernoulli_constructor(offset = y, invert = T)
val <- apply(zz, 2, function(x){
  if(all(is.na(x))) return(NA)
  func1(x) - func1(neg_grad_G) - func2(neg_grad_G)%*%(x-neg_grad_G)
})

plot(x_seq, val)

