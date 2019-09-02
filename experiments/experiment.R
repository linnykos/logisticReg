rm(list=ls())
x = 4
set.seed(x)
z <- rep(.5, 2)
f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
plane <- .plane(basis = matrix(rnorm(2), nrow = 2, ncol = 1), offset = c(-.1,.1))

# res <- .projection_bregman(z, plane, distr_class = "bernoulli")
#########

distr_class = "bernoulli"
max_iter = 100
tol = 1e-3

#res <- .setup_bregman_function(plane, distr_class)
######

stopifnot(all(dim(plane$basis) == c(2,1)))
f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
plane_list <- list(.plane(basis = matrix(c(0,1), 2, 1), offset = c(tol, tol)),
                   .plane(basis = matrix(c(1,0), 2, 1), offset = c(tol, tol)),
                   .plane(basis = matrix(c(0,1), 2, 1), offset = c(1-tol, tol)),
                   .plane(basis = matrix(c(1,0), 2, 1), offset = c(tol, 1-tol)))
point_mat <- sapply(plane_list, function(p){
  .intersect_two_lines(plane$A, p$A, plane$b, p$b)
})

idx <- which(apply(point_mat, 2, function(x){all(x >= 0) & all(x < 1)}))
stopifnot(length(idx) == 2)
