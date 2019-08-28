rm(list=ls())

func <- function(x,z){
  t(x)%*%z - sum(sapply(x, function(i){log(1+exp(i))}))
}

z = 1
x_vec <- seq(-10, 10, length.out = 100)
y_vec <- sapply(x_vec, function(x){func(x,z)})
plot(x_vec, y_vec)

plot(x_vec, log(1+exp(x_vec)))
