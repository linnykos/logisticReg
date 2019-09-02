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

res <- .setup_prox_function(plane, distr_class)
f <- res$f; grad_f <- res$grad_f; prox <- res$prox
grad_f_z <- grad_f(z)


g <- function(x){f(x) - f(z) - t(grad_f_z) %*% (x-z)}
grad_g <- function(x){grad_f(x) - grad_f_z}
G_t <- function(x, eta){(x - prox(x - eta * (grad_f(x) - grad_f_z)))/eta}

iter <- 1
x_prev <- rep(Inf, length(z))
x_current <- plane$offset + plane$basis %*% rep(1, ncol(plane$basis))

x_prev <- x_current

# eta <- .backtrack_line_search(x_current, g, grad_g, G_t)
##########

beta = 0.5
eta_init = 1

eta <- eta_init
gx <- g(x)
grad_gx <- grad_g(x)
Gtx <- G_t(x, eta)
counter <- 1



