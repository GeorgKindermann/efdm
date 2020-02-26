state0Area <- readRDS("./dat/state0Area.RData")
if(is.list(state0Area)) {
  state0Area <- "attr<-"(do.call(Matrix::sparseVector, state0Area[-1])
                      , names(state0Area)[1], state0Area[[1]])}
transition <- readRDS("./dat/transition.RData")

#Go to the next timestep
#Performance could be improved
#If needed harvests and residuals can be done for groups
retResiduals <- retHarvest <- 0
t0 <- state0Area
t1 <- t0
t1[rowSums(expand.grid(transition$allToOne$from, transition$static))] <- 0
t1[rowSums(expand.grid(unlist(lapply(transition$fromSome, "[[", 1)), transition$static))] <- 0
t1[rowSums(expand.grid(unlist(lapply(transition$fromMany$fs, "[[", 1)), transition$static))] <- 0
for(i in seq_len(nrow(transition$allToOne))) {
  j <- rowSums(expand.grid(transition$static, transition$allToOne$to[i]))
  k <- rowSums(expand.grid(transition$static, transition$allToOne$from[i]))
  t1[j] <- t1[j] + t0[k]
  retHarvest <- retHarvest + sum(t0[k] * transition$allToOne$harvest[i])
  retResiduals <- retResiduals + sum(t0[k] * transition$allToOne$residuals[i])
}
for(k in seq_along(transition$fromSome)) {
  i <- rowSums(expand.grid(transition$fromSome[[k]]$to, transition$static))
  j <- rowSums(expand.grid(transition$fromSome[[k]]$from, transition$static))
  #Currently (1.2-18) there is a bug in Matrix::sparseVector which does
  #not allow the following, so I cast it to vector what solves to problem
  #t1[i] <- t1[i] + t0[j] * transition$fromSome[[k]]$share
  t1[i] <- as.vector(t1[i]) + as.vector(t0[j]) * transition$fromSome[[k]]$share
  retHarvest <- retHarvest + sum(as.vector(t0[j]) * transition$fromSome[[k]]$share * transition$fromSome[[k]]$harvest)
  retResiduals <- retResiduals + sum(as.vector(t0[j]) * transition$fromSome[[k]]$share * transition$fromSome[[k]]$residuals)
}
i <- rowSums(expand.grid(transition$static, transition$fromMany$to))
t1[i] <- t1[i] + unlist(lapply(transition$fromMany$fs, function(j) matrix(t0[transition$static + rep(j[[1]], each=length(transition$static))], ncol=length(j[[1]])) %*% j[[2]]))
retHarvest <- retHarvest +  sum(unlist(lapply(transition$fromMany$fs, function(j) matrix(t0[transition$static + rep(j[[1]], each=length(transition$static))], ncol=length(j[[1]])) %*% (j[[2]] * j[[3]]))))
retResiduals <- retResiduals +  sum(unlist(lapply(transition$fromMany$fs, function(j) matrix(t0[transition$static + rep(j[[1]], each=length(transition$static))], ncol=length(j[[1]])) %*% (j[[2]] * j[[4]]))))
  
saveRDS(if(is(t1, "sparseVector")) {
  list(idxMul=attr(t1, "idxMul"), t1@x, t1@i, t1@length)
} else {t1}, file="./dat/t1.RData", compress="xz")

saveRDS(list(retResiduals, retHarvest), file="./dat/harvRes.RData")

