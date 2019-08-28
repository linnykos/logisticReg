set.seed(10)
dat <- matrix(rnorm(40), 4, 10)
a_idx <- c(1,4,5,6)
b_idx <- c(1:6)

######

stopifnot(all(a_idx %in% b_idx), max(b_idx) <= ncol(dat))
n <- nrow(dat); d <- ncol(dat)
if(all(b_idx %in% a_idx)) return(matrix(0, ncol = 1, nrow = n))

minusb <- c(1:d)[-b_idx]
d_minusb <- .diagonal_matrix(d, minusb)
null_d_minusb <- .nullspace(d_minusb)

d_bnota <- .diagonal_matrix(d, setdiff(b_idx, a_idx))
# compute D_{B \backslash A} (X P_{null(D_{-B})})^+
term2 <- d_bnota %*% MASS::ginv(dat %*% .projection_matrix(null_d_minusb))

null_dat <- .nullspace(dat)

# Question: how do you intersect null_d_minusb and null_dat
super_mat <- cbind(null_d_minusb, -null_dat)
null_intersection <- MASS::Null(t(super_mat))

intersection <- apply(null_intersection, 2, function(x){
  null_d_minusb %*% x[1:ncol(null_d_minusb)]
})

# check that vectors res2 are in both null_dat and null_d_minusb
tmp <- cbind(null_dat, .orthogonal_basis(null_dat))
solve(tmp, intersection[,1])

#
#
# # project each vector in null_d_minub into null_dat
# proj_mat <- null_d_minusb %*% t(null_d_minusb)
# res <- proj_mat %*% null_dat
#
# # now, determine the rank
# Matrix::rankMatrix(res)
# res2 <- .orthogonalize(res)
#

