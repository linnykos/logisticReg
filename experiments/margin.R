rm(list=ls())
X = matrix(c(1,1,2,0,2,3), ncol = 2, byrow = T)
y = c(-1,-1,1)
dat <- as.data.frame(cbind(X, y))
res <- e1071::svm(y~. , data = dat, type = "C-classification", cost = 10000000,
                  scale = FALSE, kernel = "linear", shrinking = FALSE)
res$coefs
y_alt =  as.factor(y)
lev <- levels(y_alt)
y_alt <- as.integer(y_alt)

beta = rep(0, 2)
for(i in 1:length(res$index)){
  beta <- beta + res$coefs[i]*X[res$index[i],]
}
beta_0_vec <- sapply(1:length(res$index), function(i){
  X[res$index[i],]%*%beta+y[res$index[i]]
})
range(beta_0_vec)
2/.l2norm(beta)

### hm...
t(res$coefs) %*% y[res$index]
sum(res$coefs)

## try svmpath
svmpath::svmpath(X, y)

library(e1071)
set.seed (1)
X=matrix (rnorm (150*2) , ncol =2)
X[1:100 ,]=X[1:100 ,]+10
X[101:150 ,]= X[101:150 ,] -10
y=c(rep (-1 ,100) ,rep (1 ,50) )
dat=data.frame(x=X,y=y)
res = e1071::svm(y~. , data = dat, type = "C-classification", cost = 10000000,
                scale = FALSE, kernel = "linear", shrinking = FALSE)
beta = rep(0, 2)
for(i in 1:length(res$index)){
  beta <- beta + res$coefs[i]*X[res$index[i],]
}
beta_0_vec <- sapply(1:length(res$index), function(i){
  X[res$index[i],]%*%beta+y[res$index[i]]
})
range(beta_0_vec)
2/.l2norm(beta)


length(svmfit$index)
dim(svmfit$coefs)
sum(svmfit$coefs)

plot(svmfit , dat[train ,])
