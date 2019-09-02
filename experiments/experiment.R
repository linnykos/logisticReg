rm(list=ls())
x = 1
set.seed(x)
z <- rep(.5, 2)
f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
plane <- .plane(basis = matrix(rnorm(2), nrow = 2, ncol = 1), offset = c(-.1,.1))

# res <- .projection_bregman(z, plane, distr_class = "bernoulli")
#########

distr_class = "bernoulli"
max_iter = 100
tol = 1e-3

res <- .setup_bregman_function(plane, distr_class)
