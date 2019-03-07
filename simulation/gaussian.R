rm(list=ls())
library(simulation)
library(logisticReg)

paramMat <- as.matrix(expand.grid(seq(0, 0.6, length.out = 50),
                                  seq(0, 10, length.out = 50),
                                  1000))
colnames(paramMat) <- c("kappa", "gamma", "n")
trials <- 50

################

rule <- function(vec){
  d <- max(round(vec["kappa"] * vec["n"]), 1)
  X <- MASS::mvrnorm(n = vec["n"], rep(0, d), diag(d))
  beta <- rep(sqrt(vec["gamma"])/sqrt(d), d)
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
