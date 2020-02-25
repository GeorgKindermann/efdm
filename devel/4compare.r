dat <- read.csv("./dat/inventoryData.csv.bz2")
load("./dat/t1.RData")
load("./dat/breaks.RData")
n <- nrow(dat)
c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

sum(t1) #Area t1

#Observations at t1
#t0 <- tt <- c("volume","diameter","temperature")
t0 <- tt <- c("volume","diameter")
t0[1:2]  <- paste0(t0[1:2], 1)
key <- as.data.frame(lapply(setNames(seq_along(tt),tt), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- rep("yield", 2)
t0 <- paste0(tt, c("C", "NC"))
keyY <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
stateT1 <- xtabs(area ~ ., cbind(key, yield=c(keyY$yieldC, keyY$yieldNC), species=rep(0:1, each=n), area=c(dat[,"area"] * dat$shareC1, dat[,"area"] * (1-dat$shareC1))))
#Compare for differences
summary(apply(array(t1, c(attr(t1, "idxMul")[-1], length(t1)) / attr(t1, "idxMul")), -3, sum) - as.vector(stateT1))
