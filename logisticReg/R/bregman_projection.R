.projection_euclidean <- function(point, plane){
  stopifnot(class(plane) == "plane")
  if(ncol(plane$basis) == 0) return(plane$offset)

  point2 <- point - plane$offset
  proj_mat <- .projection_matrix(plane$basis)
  proj_point <- proj_mat %*% point2

  proj_point + plane$offset
}

# apply a bregman projection onto a plane
# can only work if the plane is a line or a point
.projection_bregman <- function(z, plane, distr_class = "bernoulli", offset = rep(0,2), invert = F,
                                intersection_1 = NA, intersection_2 = NA, max_iter = 100, tol = 1e-3,
                                grid_size = 1000){
  stopifnot(ncol(plane$A) == 2, nrow(plane$A) %in% c(1,2))
  if(nrow(plane$A) == ncol(plane$A)) {
    res <- as.numeric(plane$b)
    attributes(res) <- list("iteration" = NA, "tolerance" = NA)
    return(res)
  }
  res <- .setup_bregman_function(plane, distr_class, offset = offset, invert = invert)
  point_mat <- res$point_mat
  f <- res$f; grad_f <- res$grad_f

  # check intersection points are on the plane
  if(!all(is.na(intersection_1)) | !all(is.na(intersection_2))){
    stopifnot(sum(abs(plane$A %*% intersection_1 - plane$b)) <= tol,
              sum(abs(plane$A %*% intersection_2 - plane$b)) <= tol)
  }

  # construct the grid
  if(all(is.na(intersection_1)) & all(is.na(intersection_2))){
    x_seq <- seq(0, 1, length.out = 1000)
    grid <- sapply(x_seq, function(x){point_mat[,1] + x * apply(point_mat, 1, diff)})
  } else {
    x_seq <- seq(0, 1, length.out = 1000)
    diff_vec <- c(intersection_2[1]-intersection_1[1], intersection_2[2]-intersection_1[2])
    grid <- sapply(x_seq, function(x){intersection_1 + x * diff_vec})
  }

  func1 <- .conjugate_bernoulli_constructor(offset = y, invert = T)
  func2 <- .conjugate_grad_bernoulli_constructor(offset = y, invert = T)
  val <- apply(grid, 2, function(x){
    if(all(is.na(x))) return(NA)
    func1(x) - func1(neg_grad_G) - func2(neg_grad_G)%*%(x-neg_grad_G)
  })

  as.numeric(grid[,which.min(val)])
}

.setup_bregman_function <- function(plane, distr_class, offset = rep(0,2), invert = F, tol = 1e-3){
  proj_mat <- .projection_matrix(plane$basis)

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
   point_mat <- point_mat[,idx]

  } else if(distr_class == "gaussian"){
    f <- .conjugate_gaussian_constructor(); grad_f <- .conjugate_grad_gaussian_constructor()
    x_current <- plane$offset + plane$basis %*% rep(1, ncol(plane$basis))
    point_mat <- matrix(c(-Inf, -Inf, Inf, Inf), ncol = 2)

  } else {
    stop("distr_class not properly specified")
  }

  list(f = f, grad_f = grad_f, point_mat = point_mat)
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
  counter <- 1

  while(TRUE){
    Gtx <- G_t(x, eta)
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
