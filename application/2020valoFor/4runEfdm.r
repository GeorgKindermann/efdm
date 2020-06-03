source("~/prg/efdm/efdm.r")

breaks <- readRDS("./dat/breaks.RData")
state0Area <- readRDS("./dat/state0Area.RData")
state1Area <- readRDS("./dat/state1Area.RData")
transition <- readRDS("./dat/transition.RData")
flow <- readRDS("./dat/flow.RData")


funPs <- function(tt0, title, main="", ylim=c(0,9.5e5)) {
  r0 <- as.array(drop_simple_sparse_array(rollup(tt0, c(2,3,5:9), FUN=sum)))
  r0 <- r0[,-1]
  dimnames(r0) <- list(volume=c(0,breaks$volume),speciesGroup=breaks$speciesGroup)
  colnames(r0) <- breaks$speciesGroup
  barplot(t(rowsum(r0, 100*c(0, findInterval(breaks$volume, 0:9*100)))), las=3, col=rainbow(ncol(r0)), ylim=ylim, xlab="Stock [m3/ha]", ylab="Area [ha]", main=main)
  legend("topright", rev(colnames(r0)), fill=rev(rainbow(ncol(r0))), bty="n", title=title)
}

funPf <- function(x, title="", main="", ylim=c(0,2.5e6)) {
  barplot(t(rowsum(x[-1], 6*findInterval(c(0,breaks$diameter)[x$diameter+1], 1:11*6))), col=rainbow(ncol(x)-1), ylim=ylim, xlab="Diameter [cm]", ylab="Removed Volume [m3/Year]", main=main)
  legend("topleft", rev(colnames(x[-1])), fill=rev(rainbow(ncol(x)-1)), bty="n", title=title)
}

funS <- function(x) {
  tmp <- as.array(drop_simple_sparse_array(rollup(x, c(2,3,5:9), FUN=sum)))[,-1]
  list(area=setNames(colSums(tmp), breaks$speciesGroup), stock=setNames(colSums(tmp*stock), breaks$speciesGroup))
}

stock <- c(0, filter(breaks$volume, c(.5,.5)))
stock <- c(stock[-length(stock)], breaks$volume[length(breaks$volume)])
#stock <- breaks$stock
res <- list()
t0 <- state0Area
i <- 0
repeat {
  print(i)
  pdf(paste0("/tmp/pdf",i,".pdf"))
  tt0 <- simple_sparse_array(as.matrix(efdmIndexInv(t0$i[,1], attr(t0, "dimS")))+1, t0$v, attr(t0, "dimS"))
  funPs(tt0[,,,,,,,,], paste("Year:", 2000+7*i), "All")
  funPs(tt0[,,,,2,,,,], paste("Year:", 2000+7*i), "<= 200Ha")
  funPs(tt0[,,,,3,,,,], paste("Year:", 2000+7*i), "> 200Ha")
  res[[i+1]] <- list(as=cbind(funS(tt0[,,,,2,,,,]), funS(tt0[,,,,3,,,,])))
  saveRDS(res, file=paste0("res", i, ".RData"))
  if(i==15) break
  rm(tt0)
  tt <- efdmNextArea(t0, transition)
  t0_1f <- efdmNextFlow(tt) * tt$area / 7
  x <- aggregate(t0_1f, data.frame(efdmIndexInv(tt$from, attr(tt, "dimS"))[,c("diameter","owner")]), sum)
  funPf(aggregate(.~diameter, x[-2], sum), paste0("Year: ", 2000+7*i, "-", 2000+7*(i+1)), "All", c(0,4e6))
  funPf(x[x$owner==1, -2], paste0("Year: ", 2000+7*i, "-", 2000+7*(i+1)), "<= 200Ha", c(0,4e6))
  funPf(x[x$owner==2, -2], paste0("Year: ", 2000+7*i, "-", 2000+7*(i+1)), "> 200Ha", c(0,4e6))
  dev.off()
  res[[i+1]]["f"] <- list(aggregate(t0_1f, data.frame(efdmIndexInv(tt$from, attr(tt, "dimS"))[,c("speciesGroup","owner")]), sum))
  t0 <- efdmNextState(t0, tt)
  rm(tt, t0_1f)
  i <- i+1
}

#saveRDS(res, file="res.RData")
res <- readRDS("res13.RData")

pdf("/tmp/efdm.pdf")

#Baumartenflaechen
fun <- function(x, main) {
  colnames(x) <- seq(2000, by=7, length.out = ncol(x))
  tt <- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  barplot(x, las=3, col=rainbow(nrow(x)), ylab="Area [ha]", main=main)
  legend("topright", inset=c(-0.3,0), rev(rownames(x)), fill=rev(rainbow(nrow(x))), bty="n")
  par(tt)
}
fun(t(do.call(rbind, lapply(res, function(x) x$as[[1,1]] + x$as[[1,2]]))), "All")
fun(t(do.call(rbind, lapply(res, function(x) x$as[[1,1]]))), "<=200Ha")
fun(t(do.call(rbind, lapply(res, function(x) x$as[[1,2]]))), ">200Ha")


#Stock
#barplot(t(do.call(rbind, lapply(res, function(x) x$as[[2,1]] + x$as[[2,2]]))))
#barplot(t(do.call(rbind, lapply(res, function(x) x$as[[2,1]]))))
#barplot(t(do.call(rbind, lapply(res, function(x) x$as[[2,2]]))))

#barplot(unlist(lapply(res, function(x) {sum(x$as[[2,1]] + x$as[[2,2]]) / sum(x$as[[1,1]] + x$as[[1,2]])})))
#barplot(unlist(lapply(res, function(x) {sum(x$as[[2,1]]) / sum(x$as[[1,1]])})))
#barplot(unlist(lapply(res, function(x) {sum(x$as[[2,2]]) / sum(x$as[[1,2]])})))

x <- data.frame(year=seq(2000, by=7, length.out = length(res))
,all=unlist(lapply(res, function(x) {sum(x$as[[2,1]] + x$as[[2,2]]) / sum(x$as[[1,1]] + x$as[[1,2]])}))
,"<=200Ha"=unlist(lapply(res, function(x) {sum(x$as[[2,1]]) / sum(x$as[[1,1]])}))
,">200Ha"=unlist(lapply(res, function(x) {sum(x$as[[2,2]]) / sum(x$as[[1,2]])}))
, check.names = FALSE)
plot(x$year, x[,2], type="n", ylim=range(x[,-1]), xlab="", ylab="Stock [m³/ha]")
for(i in 2:4) lines(x[,1], x[,i], col=i-1, lwd=3)
legend("topleft", colnames(x)[-1], col=1:3, lwd=3, bty="n")


#Removals
fun <- function(x, main) {
  tt <- seq(2000, by=7, length.out = ncol(x))
  colnames(x) <- paste(tt, tt+6, sep="-")
  tt <- par(mar=c(5.5, 4.1, 4.1, 6.1), xpd=TRUE)
  barplot(x, las=3, col=rainbow(nrow(x)), ylab="Removals [m³/ha/Year]", main=main)
  legend("topright", inset=c(-0.2,0), rev(rownames(x)), fill=rev(rainbow(nrow(x))), bty="n")
  par(tt)
}
fun(t(do.call(rbind, lapply(res[-length(res)], function(x) {colSums(x$f[,3:5]) / sum(x$as[[1,1]] + x$as[[1,2]])}))), "All")
fun(t(do.call(rbind, lapply(res[-length(res)], function(x) {colSums(x$f[x$f$owner==1,3:5]) / sum(x$as[[1,1]])}))), "<=200Ha")
fun(t(do.call(rbind, lapply(res[-length(res)], function(x) {colSums(x$f[x$f$owner==2,3:5]) / sum(x$as[[1,2]])}))), ">200Ha")

dev.off()















ylim <- c(0,2.5e6)
main <- "a"
title=""
barplot(t(rowsum(x[-1], 6*findInterval(c(0,breaks$diameter)[x$diameter+1], 1:11*6))), col=rainbow(ncol(x)-1), ylim=ylim, xlab="Diameter [cm]", ylab="Removed Volume [m3/Year]", main=main)
legend("topleft", rev(colnames(x[-1])), fill=rev(rainbow(ncol(x)-1)), bty="n", title=title)



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
