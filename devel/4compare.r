dat <- read.csv("./dat/inventoryData.csv.bz2")
t1 <- readRDS("./dat/t1.RData")
if(is.list(t1)) {t1 <- "attr<-"(do.call(Matrix::sparseVector, t1[-1]), names(t1)[1], t1[[1]])}
state01Stock <- readRDS("./dat/state01Stock.RData")
if(is.list(state01Stock)) {state01Stock <- "attr<-"(do.call(Matrix::sparseVector, state01Stock[-1]), names(state01Stock)[1], state01Stock[[1]])}
breaks <- readRDS("./dat/breaks.RData")
harvRes <- readRDS("./dat/harvRes.RData")
n <- nrow(dat)
c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

#Area
sum(t1)
sum(dat$area)

#Stock
sum(apply(array(t1, c(attr(t1, "idxMul")[-1], length(t1)) / attr(t1, "idxMul")), -3, sum) * apply(array(state01Stock, c(attr(state01Stock, "idxMul")[-1], length(state01Stock)) / attr(state01Stock, "idxMul")), -3, max))
sum(dat$area * (dat$volume0 + dat$volume1)/2)

#Harvest
harvRes[[2]]
sum(dat$area * dat$harvest)

#Residuals
harvRes[[1]]
sum(dat$area * dat$residuals)

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
