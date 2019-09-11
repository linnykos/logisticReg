rm(list=ls())
set.seed(100)
X = matrix(rnorm(100), ncol = 2)
y = runif(50)

fit = glmnet::glmnet(x = X, y = y, family = "binomial", intercept = F)
coef(fit, s = 1.5)
