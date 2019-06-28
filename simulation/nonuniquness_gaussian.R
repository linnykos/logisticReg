rm(list=ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 5

seq_val <- seq(-5, 5, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

dist_vec <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  .distance_point_to_set(dat, as.numeric(y_mat[i,]), lambda)
})

# form into grid
dist_mat <- matrix(dist_vec, ncol = length(seq_val))
colnames(dist_mat) <- seq_val
rownames(dist_mat) <- seq_val

image(seq_val, seq_val, dist_mat, asp = T)

bool_vec <- dist_vec <= 1e-1
bool_mat <- matrix(bool_vec, ncol = length(seq_val))

png("../figures/gaussian_nullset.png", width = 1600, height = 1600, units = "px",
    res = 300)
image(seq_val, seq_val, bool_mat, asp = T)
graphics.off()

#############################################

rm(list=ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})
lambda <- 3
seq_val <- seq(-5, 5, length.out = 100)
y_mat <- expand.grid(seq_val, seq_val)

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

png("../figures/gaussian_regions.png", width = 1600, height = 1600, units = "px",
    res = 300)
plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)
graphics.off()
