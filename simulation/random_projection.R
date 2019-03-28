rm(list=ls())
library(simulation)
library(logisticReg)

grid_size <- 30
paramMat <- as.matrix(expand.grid(seq(0, 0.95, length.out = grid_size),
                                  seq(0, 0.6, length.out = grid_size),
                                  1000))
colnames(paramMat) <- c("k_percentage", "kappa", "n")
trials <- 50

################

rule <- function(vec){
  n <- vec["n"]
  d <- max(round(vec["kappa"] * n), 1)
  sigma <- matrix(0.5, d, d)
  diag(sigma) <- 1
  X <- MASS::mvrnorm(n = n, rep(0, d), sigma)

  # do a nonparanomral transformation
  X <- apply(X, 2, function(x){
    sign(x)*abs(x)^0.5
  })

  list(X = X, y = y)
}


criterion <- function(dat, vec, y){
  n <- nrow(dat$X); d <- ncol(dat$X)
  k <- max(round(vec["k_percentage"]*d), 1)
  Z <- logisticReg::random_projection(dat$X, k)

  quant_vec <- logisticReg::gaussian_check(Z, prob_vec = c(0.5, 1))

  tmp <-  as.matrix(dist(dat$Z))/as.matrix(dist(dat$X))
  diag(tmp) <- 1
  distortion_vec <- range(tmp)

  list(quant_vec = quant_vec, distortion_vec = distortion_vec)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[2,]), paramMat[2,], 1)
# set.seed(1); criterion(rule(paramMat[900,]), paramMat[900,], 1)

#################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "../results/random_projection_tmp.RData",
                                        verbose = T)

save.image("../results/random_projection.RData")
