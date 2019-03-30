rm(list=ls())
load("../results/random_projection.RData")

ratio_vec <- sapply(res, function(x){
  tmp <- sapply(x, function(y){y$quant_vec[1]})
  median(tmp)
})

kappa_vec <- seq(0, 0.6, length.out = grid_size)
k_percentage <- seq(0, 0.95, length.out = grid_size)
ratio_mat <- t(matrix(ratio_vec, ncol = grid_size))
colnames(ratio_mat) <- k_percentage
rownames(ratio_mat) <- kappa_vec

image(x = kappa_vec, y = k_percentage,
      z = ratio_mat, breaks = seq(0, max(ratio_mat), length.out = 50),
      col = gray.colors(49),
      asp = max(kappa_vec)/max(k_percentage), xlab = "Kappa", ylab = "K/p")
# basically requires exactly gaussian data to kick in...

distortion_vec <- sapply(res, function(x){
  tmp <- sapply(x, function(y){abs(y$distortion_vec - 1)})
  median(tmp)
})

kappa_vec <- seq(0, 0.6, length.out = grid_size)
k_percentage <- seq(0, 0.95, length.out = grid_size)
ratio_mat <- t(matrix(distortion_vec, ncol = grid_size))
colnames(ratio_mat) <- k_percentage
rownames(ratio_mat) <- kappa_vec

image(x = kappa_vec, y = k_percentage,
      z = ratio_mat, breaks = seq(0, max(distortion_vec), length.out = 50),
      col = gray.colors(49),
      asp = max(kappa_vec)/max(k_percentage), xlab = "Kappa", ylab = "K/p")
