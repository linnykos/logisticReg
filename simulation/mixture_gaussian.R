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
  idx <- sample(c(1,2), n, replace = T)

  mean_vec1 <- c(-2, rep(0, d-1))
  mean_vec2 <- c(2, rep(0, d-1))

  X <- rbind(MASS::mvrnorm(n = length(which(idx == 1)), mean_vec1, diag(d)),
             MASS::mvrnorm(n = length(which(idx == 2)), mean_vec2, diag(d)))

  grand_mean <- .5*mean_vec1 + .5*mean_vec2
  pop_covariance <- (.5*diag(d) + .5*diag(d)) + (.5*(mean_vec1 - grand_mean)%*%t(mean_vec1 - grand_mean) +
                                                   .5*(mean_vec2 - grand_mean)%*%t(mean_vec2 - grand_mean))

  # load all the signal on the first coordinate
  beta <- rep(vec["gamma"]/sqrt(d), d)

  # OH ... wait. Now I see why non-Gaussian is hard...


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
                                        filepath = "../results/mixture_gaussian_tmp.RData",
                                        verbose = T)

save.image("../results/mixture_gaussian.RData")
