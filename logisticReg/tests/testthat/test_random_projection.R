context("Test random projection")

## random_projection is correct

test_that("random_projection works", {
  set.seed(10)
  n <- 100; d <- 500
  X <- MASS::mvrnorm(n, rep(0,d), diag(d))
  beta <- 1
  epsilon <- .5
  k <- round(log(n)*(4+2*beta)/(epsilon^2/2-epsilon^3/3))

  Z <- random_projection(X, k)

  expect_true(all(dim(Z) == c(n, k)))
  expect_true(is.matrix(Z))
})

################

## check_pairwise_distances is correct

test_that("check_pairwise_distances works", {
  set.seed(10)
  n <- 100; d <- 500
  X <- MASS::mvrnorm(n, rep(0,d), diag(d))
  beta <- 1
  epsilon <- .5
  k <- round(log(n)*(4+2*beta)/(epsilon^2/2-epsilon^3/3))

  Z <- random_projection(X, k)

  bool <- check_pairwise_distances(X, Z, epsilon = epsilon)

  expect_true(is.logical(bool))
  expect_true(bool)
})
