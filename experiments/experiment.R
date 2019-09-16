rm(list = ls())
set.seed(10)
dat <- matrix(rnorm(6), ncol = 3, nrow = 2)
dat <- apply(dat, 2, function(x){x/.l2norm(x)})

lambda <- 0.25

###############3

stopifnot(all(dim(dat) == c(2,3)))

n <- nrow(dat); d <- ncol(dat)

# populate the polytope
point_mat <- .populate_polytope(dat, lambda)
plot(point_mat[,1], point_mat[,2], asp = T, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))

polytope <- form_polytope(dat, lambda)

for(i in 1:length(polytope)){
  zz = rnorm(1000)*10
  val = sapply(1:length(zz), function(x){polytope[[i]]$plane$basis*zz[x]+polytope[[i]]$plane$offset})
  lines(val[1,], val[2,], col = "green", lwd = 2)
}

for(i in 1:length(polytope)){
  points(polytope[[i]]$intersection_1[1], polytope[[i]]$intersection_1[2], col = "red", pch = 16)
  points(polytope[[i]]$intersection_2[1], polytope[[i]]$intersection_2[2], col = "red", pch = 16)
}

y <- c(0.8, 0.9)
res <- bernoulli_solver(dat, y, lambda)
points(res[1]-.5, res[2]-.5, col = "blue", pch = 16)
points(y[1]-.5, y[2]-.5, col = "red", pch = 16)
