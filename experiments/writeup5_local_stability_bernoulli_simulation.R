rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

seq_val <- seq(0, 1, length.out = 100)
seq_val <- seq_val[-c(1, length(seq_val))]
y_mat <- expand.grid(seq_val, seq_val)

dist_vec5 <- sapply(1:nrow(y_mat), function(i){
  if(i %% floor(nrow(y_mat)/10) == 0) cat('*')
  .distance_to_stability_set5 (dat, as.numeric(y_mat[i,]), lambda, distr_class = "bernoulli")
})

# bool_vec <- (dist_vec5 <= 1e-1/5)
bool_vec <- (dist_vec5 <= 1e-1/3)

png("../figures/writeup5_bernoulli_lambda025_localstability.png", width = 1600, height = 1600, units = "px",
    res = 300)
plot(y_mat[,1], y_mat[,2], pch = 16, col = c("gray85", "red")[bool_vec+1], asp = T)
graphics.off()

save.image("../experiments/writeup5_bernoulli.RData")

################

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

for(i in 1:length(polytope)){
  zz = rnorm(1000)*10
  val = sapply(1:length(zz), function(x){polytope[[i]]$plane$basis*zz[x]+polytope[[i]]$plane$offset})
  lines(val[1,], val[2,], col = 1, lwd = 2)
}
