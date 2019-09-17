# the offset is when we are computing G^* instead of \psi^*
# assumes x and offset are length 2
.conjugate_bernoulli_constructor <- function(offset = rep(0, 2), invert = F, tol = 1e-4){
  func <- function(x){
    if(invert) xprime <- -x else xprime <- x
    if(any(xprime + offset < tol/2) | any(xprime + offset >= 1-tol/2)) return(Inf)

    xprime <- xprime + offset

    idx <- which(xprime > tol)
    if(length(idx) == 0) return(0)

    sum(sapply(xprime[idx], function(x_i){
      log(x_i/(1-x_i))*x_i - log(1/(1-x_i))
    }))
  }

  mat <- t(sapply(1:length(offset), function(i){
    vec <- c(0-offset[i], 1-tol-offset[i])
    if(invert) vec <- -vec
    sort(vec)
  }))
  colnames(mat) <- c("Lower", "Upper")
  attr(func, "domain") = mat

  func
}

.conjugate_grad_bernoulli_constructor <- function(offset = rep(0, 2), invert = F, tol = 1e-4){
  func <- function(x){
    if(invert) xprime <- -x else xprime <- x
    stopifnot(all(xprime + offset >= tol/2), all(xprime + offset <= 1-tol/2))

    xprime <- xprime + offset

    sapply(xprime, function(x_i){log(x_i/(1-x_i))})
  }

  mat <- t(sapply(1:length(offset), function(i){
    vec <- c(tol-offset[i], 1-tol-offset[i])
    if(invert) vec <- -vec
    sort(vec)
  }))
  colnames(mat) <- c("Lower", "Upper")
  attr(func, "domain") = mat

  func
}

.conjugate_gaussian_constructor <- function(){
  func <- function(x){.l2norm(x)^2/2}

  attr(func, "domain") = c(-Inf, Inf)
  func
}

.conjugate_grad_gaussian_constructor <- function(){
  func <- function(x){x}

  attr(func, "domain") = c(-Inf, Inf)
  func
}
