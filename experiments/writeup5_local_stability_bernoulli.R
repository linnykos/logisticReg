rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 3
y <- c(1,0)
glmnet::glmnet(x = dat, y = y, family = "binomial", intercept = F, standardize = F)
