rm(list=ls())
load("../results/gaussian.RData")

# large numbers = likely to exist
ratio_vec <- sapply(res, function(x){
  if(all(is.na(x))) return(NA)
  tmp <- sapply(x, function(y){if("bool" %in% names(y)) return(y$bool) else return(NA)})
  tmp <- tmp[!is.na(tmp)]
  sum(tmp)/length(tmp)
})

kappa_vec <- seq(0, 0.6, length.out = grid_size)
gamma_vec <- seq(0, 10, length.out = grid_size)
ratio_mat <- matrix(ratio_vec, ncol = grid_size)
colnames(ratio_mat) <- gamma_vec
rownames(ratio_mat) <- kappa_vec

vec = sapply(seq(0, 10, length.out = 100), function(x){h_mle(0, x)})

#prob of 0 = black
png("../figures/gaussian_simulation.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(x = kappa_vec, y = gamma_vec,
      z = ratio_mat, breaks = seq(0, 1, length.out = 50),
      col = gray.colors(49), #col = colorRampPalette(c("black", "white"))(49),
      asp = max(kappa_vec)/max(gamma_vec), xlab = "Kappa", ylab = "Gamma")

points(vec, seq(0, 10, length.out = 100), col = "red", pch = 16)
graphics.off()

########

# now plot the magnitude of the margin
margin_vec <- sapply(res, function(x){
  if(all(is.na(x))) return(NA)
  tmp <- sapply(x, function(y){if("margin" %in% names(y)) return(y$margin) else return(NA)})
  tmp <- tmp[!is.na(tmp)]
  median(tmp)
})
margin_mat <- matrix(margin_vec, ncol = grid_size)

## from https://stackoverflow.com/questions/13355176/gradient-legend-in-base

colfunc <- grDevices::colorRampPalette(c( "yellow",  "blue", "white", "red"))
margin_mat[is.na(margin_mat)] <- -1

png("../figures/gaussian_simulation_margin.png",
    height = 1500, width = 2000, res = 300, units = "px")
layout(matrix(1:2,ncol=2), width = c(3.5,1),height = c(1,1))
image(x = kappa_vec, y = gamma_vec,
      z = margin_mat, breaks = c(-1,seq(0, max(margin_mat), length.out = 50)),
      col = c(gray.colors(49)[49],colfunc(49)),
      asp = max(kappa_vec)/max(gamma_vec), xlab = "Kappa", ylab = "Gamma")

points(vec, seq(0, 10, length.out = 100), col = "red", pch = 16)

par(mar = c(0.5, 0.5, 4, 0.5))
legend_image <- grDevices::as.raster(matrix(rev(colfunc(49)), ncol=1))
plot(c(0,2),c(0,max(margin_mat)),type = 'n', axes = F,xlab = '', ylab = '',
     main = 'Margin value', cex.main = 0.8)
text(x=1.5, y = seq(0,max(margin_mat),l=7), labels = round(seq(0,max(margin_mat),l=7),2),
     cex = 0.8)
rasterImage(legend_image, 0, 0, 1,1)
graphics.off()
