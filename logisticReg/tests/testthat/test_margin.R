context("Test margin")

## margin is correct

test_that("margin works", {
  set.seed(10)
  n <- 1000; d <- 700
  X <- MASS::mvrnorm(n, rep(0,d), diag(d))
  beta <- rep(10/sqrt(d), d)
  y <- generate_y_from_x(X, beta = beta)
  bool <- existence(X, y) # should return FALSE, meaning there is a margin
  stopifnot(bool == FALSE)

  res <- margin(X, y)
  expect_true(is.numeric(res))
  expect_true(res > 0)
  expect_true(length(res) == 1)
})

test_that("margin gives the correct value", {
  X = matrix(c(1,1,2,0,2,3), ncol = 2, byrow = T)
  y = c(-1,-1,1)
  bool <- existence(X, y) # should return FALSE, meaning there is a margin
  stopifnot(bool == FALSE)

  res <- margin(X, y)
  expect_true(abs(res - sqrt(5)) <= 1e-6)
})

test_that("margin works for many instances",{
  trials <- 15
  res_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    n <- 1000; d <- round(.15*n)
    X <- MASS::mvrnorm(n, rep(0,d), diag(d))
    beta <- rep(10/sqrt(d), d)
    y <- generate_y_from_x(X, beta = beta)
    bool <- existence(X, y)
    if(!bool){
      margin(X, y)
    } else {
      0
    }
  })

  expect_true(all(res_vec >= 0))
})


