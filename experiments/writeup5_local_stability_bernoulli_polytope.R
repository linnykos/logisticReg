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

# shift polytope
for(i in 1:length(polytope)){
  tmp <- .plane(basis = polytope[[i]]$plane$basis,
                offset = polytope[[i]]$plane$offset + rep(0.5, 2))
  attr(tmp, "model") <- attr(polytope[[i]]$plane, "model")
  polytope[[i]]$plane <- tmp
  polytope[[i]]$intersection_1 <- polytope[[i]]$intersection_1 + rep(0.5, 2)
  polytope[[i]]$intersection_2 <- polytope[[i]]$intersection_2 + rep(0.5, 2)
}

########

model_mat <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  res <- bernoulli_solver (dat, as.numeric(y_mat[i,]), lambda, polytope = polytope)
  res$model
})

sign_string <- apply(model_mat, 2, function(x){paste0(x, collapse = "_")})
table(sign_string)
sign_string <- as.numeric(as.factor(sign_string))
plot(y_mat[,1], y_mat[,2], pch = 16, col = sign_string, asp = T)
