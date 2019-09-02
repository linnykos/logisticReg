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
    plane <- .plane(basis = matrix(rnorm(2), nrow = 2, ncol = 1), offset = c(-.1,.1))

    res <- .projection_bregman(z, plane, distr_class = "bernoulli")

    sum(abs(plane$A %*% res - plane$b)) <= 1e-6
  })



  # build a grid of points

})


test_that(".projection_bregman is minimized along a grid of points", {
  z <- rep(.5, 2)
  f <- .conjugate_bernoulli; grad_f <- .conjugate_grad_bernoulli
  plane <- .plane(basis = matrix(c(.7,.3), nrow = 2, ncol = 1), offset = c(-.1,.1))

  res <- .projection_bregman(z, f, grad_f, plane)

  # build a grid of points

})
