rm(list=ls())
library(simulation)
library(logisticReg)

grid_size <- 30
paramMat <- as.matrix(expand.grid(seq(0, 0.6, length.out = grid_size),
                                  seq(0, 10, length.out = grid_size),
                                  1000))
colnames(paramMat) <- c("kappa", "gamma", "n")
trials <- 50

################

rule <- function(vec){
  n <- vec["n"]
  d <- max(round(vec["kappa"] * n), 1)
  X <- MASS::mvrnorm(n = n, rep(0, d), diag(d))
  beta <- rep(vec["gamma"]/sqrt(d), d)
  y <- logisticReg::generate_y_from_x(X, beta_0 = 0, beta = beta)

  list(X = X, y = y)
}

criterion <- function(dat, vec, y){
  logisticReg::existence(dat$X, dat$y)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[2,]), paramMat[2,], 1)
# set.seed(1); criterion(rule(paramMat[2500,]), paramMat[2500,], 1)

#################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "../results/gaussian_tmp.RData",
                                        verbose = T)

save.image("../results/gaussian.RData")
