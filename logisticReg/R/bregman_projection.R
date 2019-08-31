.projection_euclidean <- function(point, plane){
  stopifnot(class(plane) == "plane")
  if(ncol(plane$basis) == 0) return(plane$offset)

  point2 <- point - plane$offset
  proj_mat <- .projection_matrix(plane$basis)
  proj_point <- proj_mat %*% point2

  proj_point + plane$offset
}

# apply a bregman projection onto a plane
.projection_bregman <- function(z, f, grad_f, plane, max_iter = 100, tol = 1e-3){
  proj_mat <- .projection_matrix(plane$basis)
  grad_f_z <- grad_f(z)

  prox <- function(x){proj_mat %*% x}
  g <- function(x){f(x) - f(z) - t(grad_f_z) %*% (x-z)}
  grad_g <- function(x){grad_f(x) - grad_f_z}
  G_t <- function(x, eta){(x - prox(x - eta * (grad_f(x) - grad_f_z)))/eta}

  iter <- 1
  x_prev <- rep(Inf, length(z))
  x_current <- plane$offset + plane$basis %*% rep(1, ncol(plane$basis))
  while(iter < max_iter & (is.na(tol) || .l2norm(x_prev - x_current) > tol)){
    x_prev <- x_current

    eta <- .backtrack_line_search(x_current, g, grad_g, G_t)
    x_current <- prox(x_prev - eta * grad_g(x_prev))

    iter <- iter + 1
    # print(x_current)
  }

  res <- as.numeric(x_current)
  attributes(res) <- list("iteration" = iter, "tolerance" = .l2norm(x_prev - x_current))
  res
}

.backtrack_line_search <- function(x, g, grad_g, G_t, beta = 0.5, eta_init = 1, tol = 1e-3){
  eta <- eta_init
  gx <- g(x)
  grad_gx <- grad_g(x)
  Gtx <- G_t(x, eta)
  counter <- 1

  while(TRUE){
    val1 <- g(x - eta*Gtx)
    val2 <- gx - eta * t(grad_gx)%*%Gtx + eta*.l2norm(Gtx)^2/2
    if(val2 > val1) break()
    if(abs(eta) <= tol) {eta <- 0; break()}
    # print(paste0(counter, ": ", eta, " // ", val1, " vs. ", val2))
    eta <- eta*beta

    counter <- counter + 1
  }

  eta
}
