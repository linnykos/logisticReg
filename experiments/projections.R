x = matrix(c(1,5, -1,10, 2,6,
            1,-5, 10,0, 2,-3), ncol = 2, byrow = T)
plot(x[,1], x[,2], col = c(rep("blue", 3), rep("red", 3)), asp = T, pch = 16)
lines(c(-100,100), c(-100,100), col = "red", lty = 2)

y = c(rep(1,3), rep(-1,3))
beta = c(-1, 1); beta = beta/.l2norm(beta)

stopifnot(y * (x%*%beta) >= 0)

P = beta %*% t(beta)
x2 = x %*% P

points(x2[,1], x2[,2], pch = 1, col = c(rep("blue", 3), rep("red", 3)))
