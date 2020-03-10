source("efdm.r")

dat <- read.csv("./dat/inventoryData.csv.bz2")
t1 <- readRDS("./dat/t1.RData")
state01Stock <- readRDS("./dat/state01Stock.RData")
breaks <- readRDS("./dat/breaks.RData")
flow01Harv <- readRDS("./dat/harvest01.RData")
flow01Resid <- readRDS("./dat/residuals01.RData")
n <- nrow(dat)

#Area
sum(t1)
sum(dat$area)

range(xtabs(area ~ ., with(dat, data.frame(area, efdmClassGet(dat[c("volume1", "diameter1")], FALSE, breaks)))) - apply(array(t1, attr(t1, "dimS")), 1:2, sum))

tt <- cbind(dat[,c("volume1","diameter1","temperature")], yield=c(dat$yieldC, dat$yieldNC), species=rep(0:1, each=nrow(dat)), area=c(dat[,"area"] * dat$shareC1, dat[,"area"] * (1-dat$shareC1)))
range(xtabs(area ~ ., cbind(tt["area"], efdmClassGet(tt[c("yield", "species")], FALSE, breaks))) - apply(array(t1, attr(t1, "dimS")), 4:5, sum))
range(xtabs(area ~ ., cbind(tt["area"], efdmClassGet(tt[c("volume1", "diameter1", "yield", "species")], FALSE, breaks))) - apply(array(t1, attr(t1, "dimS")), c(1:2,4:5), sum))
range(xtabs(area ~ ., cbind(tt["area"], efdmClassGet(tt[c("temperature")], FALSE, breaks))) - apply(array(t1, attr(t1, "dimS")), 3, sum))


#Stock
sum(apply(array(t1, attr(t1, "dimS")), 1, sum) * state01Stock)
sum(dat$area * dat$volume1)

apply(array(t1, attr(t1, "dimS")), 1, sum) * state01Stock
xtabs(stock ~ ., with(dat, data.frame(stock=area*(volume0+volume1)/2, efdmClassGet(dat[c("volume1")], FALSE, breaks))))

apply(array(t1, attr(t1, "dimS")), 1:2, sum) * state01Stock
xtabs(stock ~ ., with(dat, data.frame(stock=area*(volume0+volume1)/2, efdmClassGet(dat[c("volume1","diameter1")], FALSE, breaks))))

#Harvest
sum(flow01Harv)
sum(dat$area * dat$harvest)

tt <- cbind(dat[,c("volume0","diameter0","temperature")], yield=c(dat$yieldC, dat$yieldNC), species=rep(0:1, each=nrow(dat)), harvest=c(dat[,"area"] * dat$shareC0 * dat$harvest, dat[,"area"] * (1-dat$shareC0) * dat$harvest))

range(xtabs(harvest ~ ., cbind(tt["harvest"], efdmClassGet(tt[c("volume0", "diameter0", "yield", "species")], FALSE, breaks))) - apply(array(flow01Harv, attr(t1, "dimS")), c(1:2,4:5), sum))


#Residuals
sum(flow01Resid)
sum(dat$area * dat$residuals)

tt <- cbind(dat[,c("volume0","diameter0","temperature")], yield=c(dat$yieldC, dat$yieldNC), species=rep(0:1, each=nrow(dat)), residuals=c(dat[,"area"] * dat$shareC0 * dat$residuals, dat[,"area"] * (1-dat$shareC0) * dat$residuals))

range(xtabs(residuals ~ ., cbind(tt["residuals"], efdmClassGet(tt[c("volume0", "diameter0", "yield", "species")], FALSE, breaks))) - apply(array(flow01Resid, attr(t1, "dimS")), c(1:2,4:5), sum))



object.size(t1) #66048 bytes
rm(res)
gc(reset = TRUE)
res <- apply(array(t1, attr(t1, "dimS")), 1:2, sum)
gc() #15.4 Mb

library(slam)
tt <- simple_sparse_array(as.matrix(efdmIndexInv(t1@i, attr(t1, "dimS")))+1, t1@x, attr(t1, "dimS"))
object.size(tt) #114256 bytes
rm(res)
gc(reset = TRUE)
res <- drop_simple_sparse_array(rollup(tt, 3:5, FUN=sum))
gc() #0.4 Mb
as.array(res)

library(spray)
tt <- spray(as.matrix(efdmIndexInv(t1@i, attr(t1, "dimS"))), t1@x)
object.size(tt) #113432 bytes
rm(res)
gc(reset = TRUE)
res <- asum(tt, 3:5)
gc() #0.6 Mb
as.array(res, TRUE)



t2 <- readRDS("./dat/t2.RData")
sum(t1)
sum(t2)
tt <- spray(as.matrix(efdmIndexInv(t1@i, attr(t1, "dimS"))), t1@x)
round(as.array(asum(tt, 3:5), TRUE))
tt <- spray(as.matrix(efdmIndexInv(t2@i, attr(t2, "dimS"))), t2@x)
round(as.array(asum(tt, 3:5), TRUE))

flow12Harv <- readRDS("./dat/harvest12.RData")
sum(flow01Harv)
sum(flow12Harv)
tt <- spray(as.matrix(efdmIndexInv(flow01Harv@i, attr(flow01Harv, "dimS"))), flow01Harv@x)
round(as.array(asum(tt, 3:5), TRUE))
tt <- spray(as.matrix(efdmIndexInv(flow12Harv@i, attr(flow12Harv, "dimS"))), flow12Harv@x)
round(as.array(asum(tt, 3:5), TRUE))

flow12Resid <- readRDS("./dat/residuals12.RData")
sum(flow01Resid)
sum(flow12Resid)
tt <- spray(as.matrix(efdmIndexInv(flow01Resid@i, attr(flow01Resid, "dimS"))), flow01Resid@x)
round(as.array(asum(tt, 3:5), TRUE))
tt <- spray(as.matrix(efdmIndexInv(flow12Resid@i, attr(flow12Resid, "dimS"))), flow12Resid@x)
round(as.array(asum(tt, 3:5), TRUE))

