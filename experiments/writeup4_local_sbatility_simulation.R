rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 3

seq_val <- seq(-5, 5, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

####################

# compute the actual partitioning based on LARS

fit_list <- lapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')

  tryCatch({
    fit = lars::lars(x = dat, y = as.numeric(y_mat[i,]),
                     normalize = F, intercept = F)
    predict(fit, s = lambda, type = "coefficients", mode = "lambda")$coefficients
  }, error = function(e){
    rep(NA, 3)
  })
})

# convert to signs
idx <- which(sapply(fit_list, function(x){all(is.na(x))}))
if(length(idx) > 0){
  fit_list <- fit_list[-idx]
  y_mat <- y_mat[-idx,]
}

sign_list <- lapply(fit_list, function(x){
  vec <- rep(0, 3)
  vec[which(x > 1e-6)] <- 1
  vec[which(x < -1e-6)] <- -1
  vec
})

sign_string <- sapply(sign_list, function(x){paste0(x, collapse = "")})
table(sign_string)
sign_string <- as.numeric(as.factor(sign_string))

plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)

######################

# draw the supposed base partition (for null model) based on hyperplanes

par(mfrow = c(1,2))
plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)

plot(NA, xlim = c(-5,5), ylim = c(-5,5), asp = T)

for(b_idx in 1:3){
  for(s_vec in c(-1,1)){
    res <- .construct_kbs(dat, b_idx, s_vec, lambda)

    new_point1 <- res$offset + res$basis*(-100)
    new_point2 <- res$offset + res$basis*(100)

    lines(c(new_point1[1], new_point2[1]), c(new_point1[2], new_point2[2]))
  }
}


###################

# compute the local stability set for formulation 1

dist_vec1 <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  .distance_to_stability_set1 (dat, as.numeric(y_mat[i,]), lambda)
})


bool_vec1 <- (dist_vec1 <= 1e-1)

par(mfrow = c(1,2))
plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)
plot(y_mat[,1], y_mat[,2], pch = 16, col = bool_vec1, asp = T)

####################


