source("~/prg/efdm/efdm.r")

breaks <- readRDS("./dat/breaks.RData")
state0Area <- readRDS("./dat/state0Area.RData")
state1Area <- readRDS("./dat/state1Area.RData")
transition <- readRDS("./dat/transition.RData")
flow <- readRDS("./dat/flow.RData")


t0 <- efdmNextState(state0Area, efdmNextArea(state0Area, transition))
tt <- efdmNextArea(t0, transition)
t1 <- efdmNextState(t0, tt)
t0_1f <- efdmNextFlow(tt) * tt$area
colSums(t0_1f)/7

f <- transition$model[[1]]
f <- transition$model[[7]]


state1AreaE <- efdmNextArea(state0Area, transition, "share")

sum(state0Area)
sum(state1Area)
sum(state1AreaE)

t0 <- simple_sparse_array(as.matrix(efdmIndexInv(state0Area@i, attr(state0Area, "dimS")))+1, state0Area@x, attr(state0Area, "dimS"))
t1 <- simple_sparse_array(as.matrix(efdmIndexInv(state1Area@i, attr(state1Area, "dimS")))+1, state1Area@x, attr(state1Area, "dimS"))
t1e <- simple_sparse_array(as.matrix(efdmIndexInv(state1AreaE@i, attr(state1AreaE, "dimS")))+1, state1AreaE@x, attr(state1AreaE, "dimS"))

r0 <- as.array(drop_simple_sparse_array(rollup(t1, c(2,3,5:9), FUN=sum)))
r0 <- r0[,-1]
dimnames(r0) <- list(volume=c(0,breaks$volume),speciesGroup=breaks$speciesGroup)
colnames(r0) <- breaks$speciesGroup
barplot(t(r0), ylim=c(0,6e5))

barplot(t(rowsum(r0, 100*c(0, findInterval(breaks$volume, 0:9*100)))), las=3, col=rainbow(ncol(r0)), ylim=c(0,7e5), xlab="Stock [m3/ha]", ylab="Area [ha]")
legend("topright", rev(colnames(r0)), fill=rev(rainbow(ncol(r0))), bty="n", title=paste("Year:", 2000))

barplot(t(rowsum(r0, 100*c(0, findInterval(breaks$volume, 0:9*100)))), las=3, col=rainbow(ncol(r0)), ylim=c(0,5e5), xlab="Stock [m3/ha]", ylab="Area [ha]", beside=TRUE)


rm(model, envir=transition)
transition$similar <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
transition$dimWgt <- c(5, 1, 1, 1, 1, 1, 10, 7, 5)
transition$absolute <- TRUE

transition$shifter <- c(0, 0, 0, 0, 0, 0, 0, 1, 1)


t0 <- state0Area
for(i in 0:15) {
  print(i)
  tt0 <- simple_sparse_array(as.matrix(efdmIndexInv(t0$i[,1], attr(t0, "dimS")))+1, t0$v, attr(t0, "dimS"))
  r0 <- as.array(drop_simple_sparse_array(rollup(tt0, c(2,3,5:9), FUN=sum)))
  r0 <- r0[,-1]
  dimnames(r0) <- list(volume=c(0,breaks$volume),speciesGroup=breaks$speciesGroup)
  colnames(r0) <- breaks$speciesGroup
  pdf(paste0("/tmp/pdf",i,".pdf"))
  barplot(t(rowsum(r0, 100*c(0, findInterval(breaks$volume, 0:9*100)))), las=3, col=rainbow(ncol(r0)), ylim=c(0,9e5), xlab="Stock [m3/ha]", ylab="Area [ha]")
  legend("topright", rev(colnames(r0)), fill=rev(rainbow(ncol(r0))), bty="n", title=paste("Year:", 2000+7*i))
  dev.off()
  t0 <- efdmNextState(t0, efdmNextArea(t0, transition))
}

i <- rowSums(expand.grid(from, static))
j <- t0[i]
head(t0@x[match(i, t0@i)])

i <- order(t0@i)
t1 <- Matrix::sparseVector(t0@x[i], t0@i[i], length(t0))

library(slam)
i <- order(t0@i)
t1 <- simple_sparse_array(t0@i[i], t0@x[i], length(t0))
i <- rowSums(expand.grid(from, static))
j <- t1[i]


tSave <- t0
t0 <- simple_sparse_array(t0@i, t0@x, length(t0))
attr(t0, "dimS") <- attr(tSave, "dimS")
t1 <- simple_sparse_array(1, 0, dim(t0))

tt <- aggregate(t0[rowSums(expand.grid(from, static))] * transition$simple[,share], list(to=rowSums(expand.grid(to, static))), sum)

t1 <- Matrix::sparseVector(0, 1, dim(t0))


state1AreaE <- efdmNextArea(state0Area, transition, "share")
t0 <- state1AreaE
share <- "share"
maxNeighbours <- 7


t0 <- efdmNextArea(state0Area, transition, "share")
t0 <- efdmNextArea(t0, transition, "share")
share <- "share"
maxNeighbours <- 7




f <- transition$model[[8]]
predict(f, from[1,])$posterior
expand.grid(to[1,])[,names(attr(t1, "dimS"))]


findInterval(fromObsI[1]-1, transition$simple$from)
findInterval(fromObsI[1], transition$simple$from)
