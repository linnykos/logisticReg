rm(list=ls())
#see if I understand how poisson regression works
set.seed(100)
X = matrix(rnorm(10000), ncol = 2)
beta = c(0.5,1)
inner_product = X %*% beta
y = sapply(1:length(inner_product), function(x){
  rpois(1, lambda = exp(inner_product[x]))
})
plot(X[,1], y, col = rgb(0,0,0,0.1), pch = 16)

fit = stats::glm.fit(x = X, y = y, family = stats::poisson(), intercept = F)
stats::coef(fit)

###########################
# now try for glmnet
set.seed(100)
X = matrix(rnorm(10000), ncol = 2)
beta = c(0.5,1)
inner_product = X %*% beta
y = sapply(1:length(inner_product), function(x){
  max(rpois(1, lambda = exp(inner_product[x]))+rnorm(1),0)
})

fit = glmnet::glmnet(x = X, y = y, family = "poisson", intercept = F)
coef(fit, s = 1.5)

##################

# let's try the grid now
set.seed(100)
X = matrix(rnorm(6), ncol = 3)
y = c(1,5)
fit = glmnet::glmnet(x = X, y = y, family = "poisson", intercept = F)
lambda = 2
coef(fit, s = lambda)

val_length <- 50; val_max <- 10
val_grid <- expand.grid(seq(0, val_max, length.out=val_length),
                        seq(0, val_max, length.out=val_length))
fit_mat <- sapply(1:nrow(val_grid), function(x){
  if(x %% floor(nrow(val_grid)/10) == 0) cat('*')

  tryCatch({
    fit = glmnet::glmnet(x = X, y = as.numeric(val_grid[x,]), family = "poisson", intercept = F)
    as.numeric(coef(fit, s = lambda))
  }, error = function(e){
    rep(NA, 4)
  })
})

fit_mat <- fit_mat[-1,]

# convert to signs
idx <- which(apply(fit_mat, 2, function(x){all(is.na(x))}))
if(length(idx) > 0){
  val_grid <- val_grid[-idx,]
  fit_mat <- fit_mat[,-idx]
}

sign_mat <- apply(fit_mat, 2, function(x){
  vec <- rep(0, 3)
  vec[which(x > 1e-6)] <- 1
  vec[which(x < -1e-6)] <- -1
  vec
})

sign_string <- apply(sign_mat, 2, function(x){paste0(x, collapse = "")})
table(sign_string)
sign_string <- as.numeric(as.factor(sign_string))
plot(val_grid[,1], val_grid[,2], pch = 16, col = sign_string, asp = T)
