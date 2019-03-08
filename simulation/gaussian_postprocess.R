rm(list=ls())
load("../results/gaussian_tmp.RData")

paramMat <- as.matrix(expand.grid(seq(0, 0.6, length.out = 50),
                                  seq(0, 10, length.out = 50),
                                  1000))
colnames(paramMat) <- c("kappa", "gamma", "n")

# large numbers = likely to exist
ratio_vec <- sapply(res, function(x){
  if(all(is.na(x))) return(NA)
  tmp <- unlist(x)
  sum(tmp)/length(tmp)
})

kappa_vec <- seq(0, 0.6, length.out = 50)
gamma_vec <- seq(0, 10, length.out = 50)
ratio_mat <- matrix(ratio_vec, ncol = 50)
colnames(ratio_mat) <- gamma_vec
rownames(ratio_mat) <- kappa_vec

#prob of 0 = black
image(x = seq(0, 0.6, length.out = 50), y = seq(0, 10, length.out = 50),
      z = ratio_mat, breaks = seq(0, 1, length.out = 50),
      col = gray.colors(49), #col = colorRampPalette(c("black", "white"))(49),
      asp = max(kappa_vec)/max(gamma_vec), xlab = "Kappa", ylab = "Gamma")

vec = sapply(seq(0, sqrt(10), length.out = 100), function(x){h_mle(0, x)})
points(vec, seq(0, 10, length.out = 100), col = "red", pch = 16)
