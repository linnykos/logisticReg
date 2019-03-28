rm(list=ls())
load("../results/mixture_gaussian_onlyfirst.RData")
res_1 <- res
load("../results/mixture_gaussian.RData")
res_m <- res
load("../results/gaussian.RData")

diff_list <- lapply(1:length(res), function(i){
  which(sapply(1:length(res[[i]]), function(x){
    if(is.na(res[[i]][[x]]) || is.na(res_1[[i]][[x]]) || is.na(res_m[[i]][[x]])) return(NA)
    res[[i]][[x]]$bool != res_1[[i]][[x]]$bool |  res[[i]][[x]]$bool != res_m[[i]][[x]]$bool
  }))
})
