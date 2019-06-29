context("Test nonunique set")

## .construct_kbs is correct

test_that(".construct_kbs works", {
  set.seed(10)
  dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
  dat <- apply(dat, 2, function(x){x/.l2norm(x)})

  b_idx <- c(1,3)
  s_vec <- c(1,-1)
  lambda <- 3

  res <- .construct_kbs(dat, b_idx, s_vec, lambda)

  expect_true(class(res) == "plane")
  expect_true(length(res) == 4)
})

test_that(".construct_kbs forms a plane with points that solve an equation", {
  set.seed(10)
  dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
  dat <- apply(dat, 2, function(x){x/.l2norm(x)})
  lambda <- 3

  trials <- 50
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    b_idx <- sample(1:3, 1)
    s_vec <- sample(c(-1,1), 1)

    res <- .construct_kbs(dat, b_idx, s_vec, lambda)

    new_point <- res$offset + res$basis*rnorm(1)
    sum(abs(t(dat[,b_idx,drop = F]) %*% new_point - lambda*s_vec)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".construct_kbs is at where the lower-dimensional intersect", {
  set.seed(10)
  dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
  dat <- apply(dat, 2, function(x){x/.l2norm(x)})
  lambda <- 3

  combn_mat <- utils::combn(c(1:3), 2)
  for(i in 1:ncol(combn_mat)){
    b_idx <- combn_mat[,i]
    sign_combn <- .powerset(c(1,2))

    for(j in 1:length(sign_combn)){
      s_vec <- rep(-1, 2)
      s_vec[sign_combn[[j]]] <- 1

      res <- .construct_kbs(dat, b_idx, s_vec, lambda)
      intersect_point <- res$offset

      # now get the planes for each model
      res1 <- .construct_kbs(dat, b_idx[1], s_vec[1], lambda)
      res2 <- .construct_kbs(dat, b_idx[2], s_vec[2], lambda)

      dis1 <- .distance_point_to_plane(intersect_point, res1)
      dis2 <- .distance_point_to_plane(intersect_point, res2)

      expect_true(abs(dis1) + abs(dis2) <= 1e-6)
    }
  }
})

