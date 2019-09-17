.projection_euclidean <- function(point, plane){
  stopifnot(class(plane) == "plane")
  if(ncol(plane$basis) == 0) return(plane$offset)

  point2 <- point - plane$offset
  proj_mat <- .projection_matrix(plane$basis)
  proj_point <- proj_mat %*% point2

  proj_point + plane$offset
}

# apply a bregman projection onto a plane
.projection_bregman <- function(z, plane, distr_class = "bernoulli", offset = rep(0,2), invert = F,  max_iter = 100, tol = 1e-3){
  if(nrow(plane$A) == ncol(plane$A)) {
    res <- as.numeric(plane$b)
    attributes(res) <- list("iteration" = NA, "tolerance" = NA)
    return(res)
  }
  res <- .setup_bregman_function(plane, distr_class, offset = offset, invert = invert)
  x_current <- res$x_current
  if(all(is.na(x_current))) {
    res <- as.numeric(x_current)
    attributes(res) <- list("iteration" = NA, "tolerance" = NA)
    return(res)
  }

  f <- res$f; grad_f <- res$grad_f; prox <- res$prox
  grad_f_z <- grad_f(z)

  g <- function(x){f(x) - f(z) - t(grad_f_z) %*% (x-z)}
  grad_g <- function(x){grad_f(x) - grad_f_z}
  G_t <- function(x, eta){(x - prox(x - eta * (grad_f(x) - grad_f_z)))/eta}

  iter <- 1
  x_prev <- rep(Inf, length(z))

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

.setup_bregman_function <- function(plane, distr_class, offset = rep(0,2), invert = F, tol = 1e-3){
  proj_mat <- .projection_matrix(plane$basis)
  prox <- function(x){(proj_mat %*% (x - plane$offset)) + plane$offset}

  if(distr_class == "bernoulli"){
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
   if(length(idx) != 2) list(f = f, grad_f = grad_f, prox = prox, x_current = rep(NA, 2))

   x_current <- point_mat[,idx[1]]

  } else if(distr_class == "gaussian"){
    f <- .conjugate_gaussian_constructor(); grad_f <- .conjugate_grad_gaussian_constructor()
    x_current <- plane$offset + plane$basis %*% rep(1, ncol(plane$basis))

  } else {
    stop("distr_class not properly specified")
  }

  list(f = f, grad_f = grad_f, prox = prox, x_current = x_current)
}

# intersect two lines of form A1%*%y=b1 and A2%*%y=b2
.intersect_two_lines <- function(A1, A2, b1, b2){
  stopifnot(all(dim(A1) == c(1,2)), all(dim(A2) == c(1,2)), length(b1) == 1, length(b2) == 1)
  A_total <- rbind(A1, A2)
  b_total <- c(b1, b2)

  solve(A_total, b_total)
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
