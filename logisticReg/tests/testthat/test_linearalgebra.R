context("Test linear algebra")

## .projection_matrix is correct

test_that(".projection_matrix works", {
  set.seed(10)
  mat <- matrix(rnorm(30),6,5)
  res <- .projection_matrix(mat)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(6,6)))
})

test_that(".projection_matrix gives a proper projection matrix", {
  set.seed(10)
  mat <- matrix(rnorm(50),10,5)
  res <- .projection_matrix(mat)

  expect_true(sum(abs(res - t(res))) <= 1e-6) #symmetric
  expect_true(sum(abs(res %*% res - res)) <= 1e-6) #idempotent
  expect_true(Matrix::rankMatrix(res) == 5) #rank

  eigen_res <- eigen(res)
  expect_true(sum(abs(eigen_res$values[1:5] - 1)) <= 1e-6) #eigenvalue
  expect_true(sum(abs(eigen_res$values[6:10])) <= 1e-6)

  expect_true(sum(abs(res %*% mat - mat)) <= 1e-6) #preserves column space
})

test_that(".projection_matrix preserves column space", {
  set.seed(20)
  mat <- matrix(rnorm(50),10,5)
  res <- .projection_matrix(mat)

  trials <- 50
  bool_vec <- sapply(1:trials, function(i){
    new_mat <- mat %*% diag(rnorm(5))
    sum(abs(res %*% new_mat - new_mat)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

########################

## .orthogonalize is correct

test_that(".orthogonalize works", {
  mat <- matrix(1:60, 10, 6)
  res <- .orthogonalize(mat)

  expect_true(all(dim(res) == c(10, 6)))
})

test_that(".orthogonalize returns normed vectors by convention", {
  set.seed(10)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonalize(mat)

  expect_true(all(apply(res, 2, .l2norm)) == 1)
})

test_that(".orthogonalize returns orthogonal vectors", {
  mat <- matrix(1:60, 10, 6)
  res <- .orthogonalize(mat)

  expect_true(sum(all(t(res)%*%res - diag(6))) <= 1e-6)
})

############################

## .orthogonal_basis is correct

test_that(".orthogonal_basis works", {
  set.seed(10)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonal_basis(mat)

  expect_true(all(dim(res) == c(10, 4)))
})

test_that(".orthogonal_basis returns orthogonal vectors", {
  set.seed(20)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonal_basis(mat)

  bool_vec <- sapply(1:ncol(res), function(i){
    sum(abs(t(mat) %*% res[,i])) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".orthogonal_basis returns normed vectors by convention", {
  set.seed(10)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonal_basis(mat)

  expect_true(all(apply(res, 2, .l2norm)) == 1)
})

test_that(".orthogonal_basis returns orthogonal vectors", {
  set.seed(30)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonal_basis(mat)

  expect_true(sum(all(t(res)%*%res - diag(4))) <= 1e-6)
})



