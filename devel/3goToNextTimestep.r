load("./dat/initstate.RData")
load("./dat/transition.RData")

#Go to the next timestep
#Performance could be improved
t0 <- stateInit
t1 <- t0
t1[rowSums(expand.grid(transition$allToOne$from, transition$static))] <- 0
t1[rowSums(expand.grid(unlist(lapply(transition$toSome, "[[", 1)), transition$static))] <- 0
t1[rowSums(expand.grid(unlist(lapply(transition$toMany$fs, "[[", 1)), transition$static))] <- 0
for(i in seq_len(nrow(transition$allToOne))) {
  j <- rowSums(expand.grid(transition$static, transition$allToOne$to[i]))
  t1[j] <- t1[j] + t0[rowSums(expand.grid(transition$static, transition$allToOne$from[i]))]
}
for(k in seq_along(transition$toSome)) {
  i <- rowSums(expand.grid(transition$toSome[[k]]$to, transition$static))
  j <- rowSums(expand.grid(transition$toSome[[k]]$from, transition$static))
  #Currently (1.2-18) there is a bug in Matrix::sparseVector which does
  #not allow the following, so I cast it to vector what solves to problem
  #t1[i] <- t1[i] + t0[j] * transition$toSome[[k]]$share
  t1[i] <- as.vector(t1[i]) + as.vector(t0[j]) * transition$toSome[[k]]$share
}
i <- rowSums(expand.grid(transition$static, transition$toMany$to))
t1[i] <- t1[i] + unlist(lapply(transition$toMany$fs, function(j) matrix(t0[transition$static + rep(j[[1]], each=length(transition$static))], ncol=length(j[[1]])) %*% j[[2]]))

save(t1, file="./dat/t1.RData", compress="xz")
