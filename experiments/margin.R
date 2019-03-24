rm(list=ls())
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

set.seed(2)
vec <- paramMat[84,]
dat <- rule(vec)

######

X = dat$X
y = dat$y
check = TRUE
tol = 1e-2

res <- e1071::svm(x = X, y = y, type = "C-classification", cost = 1000000000,
                  scale = FALSE, kernel = "linear", shrinking = FALSE)
if(check) stopifnot(all(y %in% c(-1,1)))

beta = rep(0, 2)
for(i in 1:length(res$index)){
  beta <- beta + res$coefs[i]*X[res$index[i],]
}

stopifnot(abs(sum(res$coefs)) < tol,
          length(res$coefs) == length(res$index))

#checking the distance from margin is constant
# it's unclear to me if the labels are sometimes flipped...
beta_0_vec1 <- sapply(1:length(res$index), function(i){
  X[res$index[i],]%*%beta+y[res$index[i]]
})
beta_0_vec2 <- sapply(1:length(res$index), function(i){
  X[res$index[i],]%*%beta-y[res$index[i]]
})

