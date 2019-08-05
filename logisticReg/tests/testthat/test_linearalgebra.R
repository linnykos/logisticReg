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

  expect_true(sum(abs(apply(res, 2, .l2norm) - 1)) <= 1e-6)
})

test_that(".orthogonalize returns orthogonal vectors", {
  set.seed(10)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonalize(mat)

  expect_true(sum(abs(t(res)%*%res - diag(6))) <= 1e-6)
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

  expect_true(sum(abs(apply(res, 2, .l2norm)) - 1) <= 1e-6)
})

test_that(".orthogonal_basis returns orthogonal vectors", {
  set.seed(30)
  mat <- matrix(rnorm(60), 10, 6)
  res <- .orthogonal_basis(mat)

  expect_true(sum(abs(t(res)%*%res - diag(4))) <= 1e-6)
})

##################################

## .nullspace is correct

test_that(".nullspace works", {
  set.seed(10)
  mat <- matrix(rnorm(60), 6, 10)
  res <- .nullspace(mat)

  expect_true(all(dim(res) == c(10,4)))
})

test_that(".nullspace returns normed vectors by convention", {
  set.seed(30)
  mat <- matrix(rnorm(60), 6, 10)
  res <- .nullspace(mat)

  expect_true(sum(abs(apply(res, 2, .l2norm)) - 1) <= 1e-6)
})

test_that(".nullspace returns orthogonal vectors", {
  set.seed(50)
  mat <- matrix(rnorm(60), 6, 10)
  res <- .nullspace(mat)

  expect_true(sum(abs(t(res)%*%res - diag(4))) <= 1e-6)
})

####################

## .plane is correct

test_that(".plane works", {
  set.seed(10)

  basis <- matrix(rnorm(60), 10, 6)
  offset <- rep(0, 10)
  res <- .plane(basis, offset)

  expect_true(class(res) == "plane")
  expect_true(all(sort(names(res)) == sort(c("basis", "offset", "A", "b"))))
  expect_true(all(dim(res$basis) == c(10,6)))
  expect_true(length(res$offset) == 10)
  expect_true(all(dim(res$A) == c(4,10)))
  expect_true(length(res$b) == 4)
})

# generate a point via (basis*alpha + offset) representation, and make sure
#  it satisfies Ax =b
test_that(".plane has a correct representation from basis to A", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)

    new_vec <- as.numeric(basis %*% rnorm(6)) + offset
    sum(abs(plane$A %*% new_vec - plane$b)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

# generate a point that is guaranteed to satisfy (Ax=b) called "new_vec", and make sure
#  it satisfies (basis*alpha+offset) representation by (new_vec - offset) is the basis
#  since it does not move when you project it onto the basis
test_that(".plane has a correct representation from A to basis", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)

    orth_basis <- .orthogonal_basis(t(plane$A))
    new_vec <- MASS::ginv(plane$A) %*% plane$b + orth_basis %*% rnorm(6)

    proj_mat <- .projection_matrix(plane$basis)
    new_vec_offset <- new_vec - plane$offset
    sum(abs(plane$A %*% new_vec - plane$b)) <= 1e-6 & sum(abs(proj_mat %*% new_vec_offset - new_vec_offset)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".plane agree on representation of basis and offset", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)
    plane2 <- .plane(basis = NA, offset = NA, A = plane$A, b = plane$b)

    point <- .point_on_plane(plane$A, plane$b)

    residual1 <- point - plane$offset
    basis_complete1 <- cbind(plane$basis, .orthogonal_basis(plane$basis))
    sol1 <- solve(basis_complete1, residual1)

    residual2 <- point - plane2$offset
    basis_complete2 <- cbind(plane2$basis, .orthogonal_basis(plane2$basis))
    sol2 <- solve(basis_complete2, residual2)

    all(c(dim(plane$basis) == dim(plane2$basis), all(abs(sol1[7:10]) <= 1e-6),
          all(abs(sol2[7:10]) <= 1e-6)))
  })

  expect_true(all(bool_vec))
})

test_that(".plane agree on representation of A and b", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(i){
    set.seed(10*i)
    A <- matrix(rnorm(12), 2, 6)
    b <- rnorm(2)
    plane <- .plane(basis = NA, offset = NA, A = A, b = b)
    plane2 <- .plane(basis = plane$basis, offset = plane$offset)

    point <- plane$basis %*% rnorm(4) + plane$offset

    residual1 <- plane$A %*% point - plane$b
    residual2 <- plane2$A %*% point - plane2$b

    all(abs(c(residual1, residual2)) <= 1e-6) & all(dim(plane$A) == dim(plane2$A))
  })

  expect_true(all(bool_vec))
})

#######################

test_that(".point_on_plane works", {
  set.seed(1)
  basis <- matrix(rnorm(60), 10, 6)
  offset <- rnorm(10)
  plane <- .plane(basis, offset)
  res <- .point_on_plane(plane$A, plane$b)

  expect_true(length(res) == 10)
})

test_that(".point_on_plane gives a proper point", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)

    res <- .point_on_plane(plane$A, plane$b)

    sum(abs(plane$A %*% res- plane$b)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

##########################

## .distance_point_to_plane is correct

test_that(".distance_point_to_plane works", {
  set.seed(5)
  basis <- matrix(rnorm(60), 10, 6)
  offset <- rnorm(10)
  plane <- .plane(basis, offset)
  point <- rnorm(10)
  res <- .distance_point_to_plane(point, plane)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
})

test_that(".distance_point_to_plane is always non-negative", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)
    point <- rnorm(10)

    res <- .distance_point_to_plane(point, plane)

    ifelse(res >= 0, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".distance_point_to_plane satisifies triangle inequality", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)

    point1 <- rnorm(10)
    point2 <- rnorm(10)

    res1 <- .distance_point_to_plane(point1, plane)
    res2 <- .distance_point_to_plane(point2, plane)
    dist <- .l2norm(point1 - point2)

    ifelse(abs(res2-res1) <= dist, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".distance_point_to_plane is 0 on the plane", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    basis <- matrix(rnorm(60), 10, 6)
    offset <- rnorm(10)
    plane <- .plane(basis, offset)

    point <- as.numeric(basis %*% rnorm(6)) + offset

    res <- .distance_point_to_plane(point, plane)

    abs(res) <= 1e-6
  })

  expect_true(all(bool_vec))
})

# example from https://mathinsight.org/distance_point_plane_examples
test_that(".distance_point_to_plane is correct on an example", {
  plane <- .plane(basis = NA, offset = NA, A = matrix(c(2,-2,5), nrow = 1, ncol = 3),
                  b = -8)
  point <- c(4, -4, 3)

  res <- .distance_point_to_plane(point, plane)

  expect_true(abs(res - 39/sqrt(33)) <= 1e-6)
})

# example from https://www.math.ucla.edu/~ronmiech/Calculus_Problems/32A/chap11/section5/718d63/718_63.html
test_that(".distance_point_to_plane is correct on an example 2", {
  plane <- .plane(basis = NA, offset = NA, A = matrix(c(1,-2,-2), nrow = 1, ncol = 3),
                  b = 1)
  point <- c(2, 8, 5)

  res <- .distance_point_to_plane(point, plane)

  expect_true(abs(res - 25/3) <= 1e-6)
})
