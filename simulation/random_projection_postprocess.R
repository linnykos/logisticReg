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

beta <- .5
epsilon <- .5
magic_kappa <- log(paramMat[1,"n"])*(4+2*beta)/(epsilon^2/2 -epsilon^3/3)
x_seq <- seq(0, 0.6, length.out = 1000)
y_seq <- sapply(x_seq, function(x){
  tmp <- seq(0, 0.95, length.out = 1000)
  idx_vec <- which(round(tmp*round(x*paramMat[1,"n"])) > magic_kappa)
  if(length(idx_vec) == 0) return(NA)
  min(tmp[idx_vec])
})
x_seq <- x_seq[!is.na(y_seq)]; y_seq <- y_seq[!is.na(y_seq)]


colfunc <- grDevices::colorRampPalette(c( "yellow",  "blue", "white", "red"))

png("../figures/random_projection_gaussian.png",
    height = 1500, width = 2000, res = 300, units = "px")
layout(matrix(1:2,ncol=2), width = c(3.5,1),height = c(1,1))
image(x = kappa_vec, y = k_percentage,
      z = ratio_mat, breaks = seq(0, max(ratio_mat), length.out = 50),
      col = colfunc(49),
      asp = max(kappa_vec)/max(k_percentage), xlab = "Kappa", ylab = "k/d",
      main = "Median marginal Gaussian difference")

points(x_seq, y_seq, col = "black", pch = 16)

par(mar = c(0.5, 0.5, 4, 0.5))
legend_image <- grDevices::as.raster(matrix(rev(colfunc(49)), ncol=1))
plot(c(0,2),c(min(ratio_mat),max(ratio_mat)),type = 'n', axes = F,xlab = '', ylab = '',
     main = 'Gaussian difference', cex.main = 0.8)
text(x=1.5, y = seq(min(ratio_mat),max(ratio_mat),l=7),
     labels = round(seq(min(ratio_mat),max(ratio_mat),l=7),2),
     cex = 0.8)
rasterImage(legend_image, 0, min(ratio_mat), 1,max(ratio_mat))
graphics.off()


####################

distortion_vec <- sapply(res, function(x){
  tmp <- sapply(x, function(y){abs(y$distortion_vec - 1)})
  median(tmp)
})

kappa_vec <- seq(0, 0.6, length.out = grid_size)
k_percentage <- seq(0, 0.95, length.out = grid_size)
ratio_mat <- t(matrix(distortion_vec, ncol = grid_size))
colnames(ratio_mat) <- k_percentage
rownames(ratio_mat) <- kappa_vec

beta <- .5
epsilon <- .5
magic_kappa <- log(paramMat[1,"n"])*(4+2*beta)/(epsilon^2/2 -epsilon^3/3)
x_seq <- seq(0, 0.6, length.out = 1000)
y_seq <- sapply(x_seq, function(x){
  tmp <- seq(0, 0.95, length.out = 1000)
  idx_vec <- which(round(tmp*round(x*paramMat[1,"n"])) > magic_kappa)
  if(length(idx_vec) == 0) return(NA)
  min(tmp[idx_vec])
})
x_seq <- x_seq[!is.na(y_seq)]; y_seq <- y_seq[!is.na(y_seq)]

colfunc <- grDevices::colorRampPalette(c( "yellow",  "blue", "white", "red"))


png("../figures/random_projection_distortion.png",
    height = 1500, width = 2000, res = 300, units = "px")
layout(matrix(1:2,ncol=2), width = c(3.5,1),height = c(1,1))
image(x = kappa_vec, y = k_percentage,
      z = ratio_mat, breaks = seq(0, max(distortion_vec), length.out = 50),
      col = colfunc(49),
      asp = max(kappa_vec)/max(k_percentage), xlab = "Kappa", ylab = "k/d",
      main = "Pairwise distance distortion")


points(x_seq, y_seq, col = "black", pch = 16)

par(mar = c(0.5, 0.5, 4, 0.5))
legend_image <- grDevices::as.raster(matrix(rev(colfunc(49)), ncol=1))
plot(c(0,2),
     c(min(ratio_mat),max(ratio_mat)),
     type = 'n', axes = F,xlab = '', ylab = '',
     main = 'Distortion', cex.main = 0.8)
text(x=1.5, y = seq(min(ratio_mat),max(ratio_mat),l=7),
     labels = round(seq(min(ratio_mat),max(ratio_mat),l=7),2),
     cex = 0.8)
rasterImage(legend_image, 0, min(ratio_mat), 1,max(ratio_mat))
graphics.off()
