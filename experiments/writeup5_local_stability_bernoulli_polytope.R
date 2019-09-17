#the goal is to exactly determine which points fall in which parts of the polyhedron
# in the spirit of Figure 2 of http://www.stat.cmu.edu/~ryantibs/papers/genlasuni.pdf

rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

seq_val <- seq(0, 1, length.out = 100)
seq_val <- seq_val[-c(1, length(seq_val))]
y_mat <- expand.grid(seq_val, seq_val)

polytope <- .form_polytope(dat, lambda)

########

model_mat <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  res <- bernoulli_solver (dat, as.numeric(y_mat[i,]), lambda, polytope = polytope)
  res$model
})

save.image("../experiments/writeup5_bernoulli_polytope2.RData")

sign_string <- apply(model_mat, 2, function(x){paste0(x, collapse = "_")})
table(sign_string)
sign_string <- as.numeric(as.factor(sign_string))
png("../figures/writeup5_bernoulli_lambda025_modelselection2.png", width = 1600, height = 1600, units = "px",
    res = 300)
plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)
graphics.off()

# for(i in 1:length(polytope)){
#   zz = rnorm(1000)*10
#   val = sapply(1:length(zz), function(x){polytope[[i]]$plane$basis*zz[x]+polytope[[i]]$plane$offset})
#   lines(val[1,], val[2,], col = i, lwd = 2)
# }
#
