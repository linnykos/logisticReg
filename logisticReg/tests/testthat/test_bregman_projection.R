context("Test bregman projections")

## .projection_bregman is correct

test_that(".projection_bregman works for gaussian", {
  z <- rep(.5, 2)
  f <- .conjugate_gaussian; grad_f <- .conjugate_grad_gaussian
  plane <- .plane(basis = matrix(c(.7,.3), nrow = 2, ncol = 1))

  res <- .projection_bregman(z, plane, class = "gaussian")

  expect_true(is.numeric(res))
})

test_that(".projection_bregman works for bernoulli", {
  z <- rep(.5, 2)
  f <- .conjugate_gaussian; grad_f <- .conjugate_grad_gaussian
  plane <- .plane(basis = matrix(c(.7,.3), nrow = 2, ncol = 1))

  res <- .projection_bregman(z, plane, class = "bernoulli")

  expect_true(is.numeric(res))
})

test_that(".projection_bregman results in a point on the plane", {
  trials <- 20

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    z <- rep(.5, 2)
    f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
    plane <- .plane(basis = matrix(rnorm(2), nrow = 2, ncol = 1), offset = c(.1,.1))

    res <- .projection_bregman(z, plane, distr_class = "bernoulli")
    sum(abs(plane$A %*% res - plane$b)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".projection_bregman is minimized along a grid of points", {
  trials <- 20

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    z <- rep(.5, 2)
    f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
    plane <- .plane(basis = matrix(rnorm(2), nrow = 2, ncol = 1), offset = c(.1,.1))

    res <- .projection_bregman(z, plane, distr_class = "bernoulli")
    obj_val <- .conjugate_bernoulli(res) - .conjugate_bernoulli(z) - t(grad_f(z)) %*% (res - z)

    # build a grid of points
    tol <- 1e-3
    plane_list <- list(.plane(basis = matrix(c(0,1), 2, 1), offset = c(tol, tol)),
                       .plane(basis = matrix(c(1,0), 2, 1), offset = c(tol, tol)),
                       .plane(basis = matrix(c(0,1), 2, 1), offset = c(1-tol, tol)),
                       .plane(basis = matrix(c(1,0), 2, 1), offset = c(tol, 1-tol)))
    point_mat <- sapply(plane_list, function(p){
      .intersect_two_lines(plane$A, p$A, plane$b, p$b)
    })
    idx <- which(apply(point_mat, 2, function(x){all(x > 0) & all(x < 1)}))

    point_seq <- sapply(seq(0,1,length.out=100), function(x){
      x*point_mat[,idx[1]] + (1-x)*point_mat[,idx[2]]
    })
    obj_seq <- apply(point_seq, 2, function(x){
      .conjugate_bernoulli(x) - .conjugate_bernoulli(z) - t(grad_f(z)) %*% (x - z)
    })

    min(obj_seq) - obj_val >= -1e-3
  })

  expect_true(all(bool_vec))
})
