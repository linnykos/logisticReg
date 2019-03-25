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
  d <- max(round(vec["kappa"] * n), 2)
  idx <- sample(c(1,2), n, replace = T)

  mean_vec1 <- c(-2, rep(0, d-1)); mean_vec2 <- c(2, rep(0, d-1))

  X <- rbind(MASS::mvrnorm(n = length(which(idx == 1)), mean_vec1, diag(d)),
             MASS::mvrnorm(n = length(which(idx == 2)), mean_vec2, diag(d)))

  grand_mean <- .5*mean_vec1 + .5*mean_vec2
  pop_covariance <- (.5*diag(d) + .5*diag(d)) + (.5*(mean_vec1 - grand_mean)%*%t(mean_vec1 - grand_mean) +
                                                   .5*(mean_vec2 - grand_mean)%*%t(mean_vec2 - grand_mean))

  # load all the signal on all but the first coordinate
  # beta <- c(0, sqrt(2)*vec["gamma"]/sqrt(d), rep(vec["gamma"]/sqrt(d), d-2))
  # beta <- c(vec["gamma"]/(sqrt(d)*sqrt(pop_covariance[1,1])), rep(vec["gamma"]/sqrt(d), d-1))
  beta <- c(vec["gamma"]/sqrt(pop_covariance[1,1]), rep(0, d-1))

  #stopifnot(abs(t(beta)%*%pop_covariance%*%beta - vec["gamma"]^2) < 1e-6)

  y <- logisticReg::generate_y_from_x(X, beta_0 = 0, beta = beta)

  list(X = X, y = y)
}

criterion <- function(dat, vec, y){
  bool <- logisticReg::existence(dat$X, dat$y)
  # if(!bool){
  #   margin <- logisticReg::margin(dat$X, dat$y)
  # } else {
  #   margin <- NA
  # }

  list(bool = bool) #, margin = margin)
}

# set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)
# set.seed(1); criterion(rule(paramMat[2,]), paramMat[2,], 1)
# set.seed(1); criterion(rule(paramMat[900,]), paramMat[900,], 1)

#################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = 20, as_list = T,
                                        filepath = "../results/mixture_gaussian_tmp.RData",
                                        verbose = T)

save.image("../results/mixture_gaussian_onlyfirst.RData")


## res <- res[which(sapply(res, length) > 0)]; zz <- sapply(res, function(i){length(which(sapply(i, length) == 1))}); names(zz) <- NULL; zz

