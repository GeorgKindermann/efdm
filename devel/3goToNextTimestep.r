source("efdm.r")

state0Area <- readRDS("./dat/state0Area.RData")
transition <- readRDS("./dat/transition.RData")

t0 <- state0Area
t1 <- t0
if(is.vector(state0Area)) {
  retResiduals <- retHarvest <- numeric(length(state0Area))
} else {
  retResiduals <- retHarvest <- Matrix::sparseVector(0, 1, length(state0Area))
}

static <- efdmIndexStatic(setdiff(names(attr(state0Area, "dimS")), names(transition$dimS)), attr(state0Area, "dimS"))

t1[rowSums(expand.grid(efdmIndexGet(efdmIndexInv(transition$allToOne$from, transition$dimS), attr(state0Area, "dimS")), static ))]  <- 0
t1[rowSums(expand.grid(efdmIndexGet(efdmIndexInv(unlist(lapply(transition$fromSome, "[[", 1)), transition$dimS), attr(state0Area, "dimS")), static))]  <- 0
t1[rowSums(expand.grid(efdmIndexGet(efdmIndexInv(unlist(transition$fromMany$from), transition$dimS), attr(state0Area, "dimS")), static))]  <- 0

for(i in seq_len(nrow(transition$allToOne))) {
  j <- rowSums(expand.grid(static, efdmIndexGet(efdmIndexInv(transition$allToOne$to[i], transition$dimS), attr(state0Area, "dimS"))))
  k <- rowSums(expand.grid(static, efdmIndexGet(efdmIndexInv(transition$allToOne$from[i], transition$dimS), attr(state0Area, "dimS"))))
  t1[j] <- as.vector(t1[j]) + as.vector(t0[k])
  retHarvest[k] <- as.vector(retHarvest[k]) + as.vector(t0[k]) * transition$allToOne$harvest[i]
  retResiduals[k] <- as.vector(retResiduals[k]) + as.vector(t0[k]) * transition$allToOne$residuals[i]
}

for(k in seq_along(transition$fromSome)) {
  i <- rowSums(expand.grid(efdmIndexGet(efdmIndexInv(transition$fromSome[[k]]$to, transition$dimS), attr(state0Area, "dimS")), static))
  j <- rowSums(expand.grid(efdmIndexGet(efdmIndexInv(transition$fromSome[[k]]$from, transition$dimS), attr(state0Area, "dimS")), static))
  #Currently (1.2-19) there is a bug in Matrix::sparseVector which does
  #not allow the following, so I cast it to vector what solves to problem
  #t1[i] <- t1[i] + t0[j] * transition$fromSome[[k]]$share
  t1[i] <- as.vector(t1[i]) + as.vector(t0[j]) * transition$fromSome[[k]]$share
  retHarvest[j] <- as.vector(retHarvest[j]) + as.vector(t0[j]) * transition$fromSome[[k]]$share * transition$fromSome[[k]]$harvest
  retResiduals[j] <- as.vector(retResiduals[j]) + as.vector(t0[j]) * transition$fromSome[[k]]$share * transition$fromSome[[k]]$residuals
}

i <- rowSums(expand.grid(static, efdmIndexGet(efdmIndexInv(unlist(transition$fromMany$to), transition$dimS), attr(state0Area, "dimS"))))
t1[i] <- as.vector(t1[i]) + unlist(lapply(split(transition$fromMany, seq_len(nrow(transition$fromMany))), function(j) matrix(as.vector(t0[static + rep(efdmIndexGet(efdmIndexInv(unlist(j$from), transition$dimS), attr(state0Area, "dimS")), each=length(static))]), ncol=length(unlist(j$from)))  %*% unlist(j$share)))
for(i in seq_len(nrow(transition$fromMany))) { #This is a slow Solution
  j <- rowSums(expand.grid(static, efdmIndexGet(efdmIndexInv(transition$fromMany$from[[i]], transition$dimS), attr(state0Area, "dimS"))))
  retHarvest[j] <- as.vector(retHarvest[j]) + as.vector(t0[j]) * transition$fromMany$share[[i]] * transition$fromMany$harvest[[i]]
  retResiduals[j] <- as.vector(retResiduals[j]) + as.vector(t0[j]) * transition$fromMany$share[[i]] * transition$fromMany$residuals[[i]]
}

saveRDS(t1, file="./dat/t1.RData", compress="xz")
saveRDS(retHarvest, file="./dat/harvest01.RData", compress="xz")
saveRDS(retResiduals, file="./dat/residuals01.RData", compress="xz")

