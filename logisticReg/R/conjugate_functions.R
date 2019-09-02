.conjugate_bernoulli <- function(x, tol = 1e-4){
  if(any(x < 0) | any(x >= 1-tol)) return(Inf)

  idx <- which(x > tol)
  if(length(idx) == 0) return(0)

  sum(sapply(x[idx], function(x_i){
    log(x_i/(1-x_i))*x_i - log(1/(1-x_i))
  }))
}

.conjugate_grad_bernoulli <- function(x, tol = 1e-4){
  stopifnot(all(x >= 0), all(x <= 1-tol))

  sapply(x, function(x_i){log(x_i/(1-x_i))})
}

.conjugate_gaussian <- function(x){
  .l2norm(x)^2/2
}

.conjugate_grad_gaussian <- function(x){
  x
}
