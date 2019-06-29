rm(list=ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})
lambda <- 3

plot(NA, xlim = c(-5,5), ylim = c(-5,5), asp = T)

for(b_idx in 1:3){
  for(s_vec in c(-1,1)){
    res <- .construct_kbs(dat, b_idx, s_vec, lambda)

    new_point1 <- res$offset + res$basis*(-100)
    new_point2 <- res$offset + res$basis*(100)

    lines(c(new_point1[1], new_point2[1]), c(new_point1[2], new_point2[2]))
  }
}

combn_mat <- utils::combn(c(1:3), 2)
for(i in 1:ncol(combn_mat)){
  b_idx <- combn_mat[,i]
  sign_combn <- .powerset(c(1,2))
  for(j in 1:length(sign_combn)){
    s_vec <- rep(-1, 2)
    s_vec[sign_combn[[j]]] <- 1

    res <- .construct_kbs(dat, b_idx, s_vec, lambda)
    points(res$offset[1], res$offset[2], pch = 16, col = "red")
  }
}

# plot all 3-variable models
sign_combn <- .powerset(c(1:3))
for(j in 1:length(sign_combn)){
  s_vec <- rep(-1, 3)
  s_vec[sign_combn[[j]]] <- 1

  res <- .construct_kbs(dat, 1:3, s_vec, lambda)
  print(res$offset)
  points(res$offset[1], res$offset[2], pch = 16, col = "green")
}
