rm(list=ls())
load("../results/mixture_gaussian_withfirst.RData")

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

vec <- sapply(seq(0, 10, length.out = 100), function(x){h_mle(0, x)})

#prob of 0 = black
png("../figures/mixture_gaussian_simulation_withfirst.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(x = kappa_vec, y = gamma_vec,
      z = ratio_mat, breaks = seq(0, 1, length.out = 50),
      col = gray.colors(49), #col = colorRampPalette(c("black", "white"))(49),
      asp = max(kappa_vec)/max(gamma_vec), xlab = "Kappa", ylab = "Gamma")

points(vec, seq(0, 10, length.out = 100), col = "red", pch = 16)
graphics.off()
