c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

x <- read.csv("./EFDMdat/inventoryData.csv.bz2")

y <- read.table("./res/rawoutput.txt")

tt <- aggregate(sapply(0:2, function(i) {rowSums(y[,grep(paste0("step",i,"."), colnames(y))])}), y[,c("volume", "diameter")], FUN=sum)

plot(tt[,-2:-1], pch=".")

cbind(tt[,1:2], round(tt[,-2:-1]/sum(tt[,3])*100,2))





#t0 <- with(x, aggregate(list(area=c(area*antc, area*(1-antc))), data.frame(volume=c(v0ck,v0nck), diameter=c(d0ck, d0nck), siteIndex=sik, species=rep(c("coniferous", "nonConiferous"), each=length(area)), slope=slpk, transportDistance=distk, owner=ownerk, harvestType=harvk), FUN=sum))

o0 <- with(x, aggregate(list(areaO=c(area*antc, area*(1-antc))), data.frame(volume=c(v0ck,v0nck), diameter=c(d0ck, d0nck), siteIndex=sik, species=rep(c("coniferous", "nonConiferous"), each=length(area)), slope=slpk, transportDistance=distk, owner=ownerk), FUN=sum))

e0 <- cbind(y[,grep("^step[0-9]+\\.", names(y), invert=TRUE)], areaE=rowSums(y[,grep("step0\\.", colnames(y))]))

me <- merge(o0, e0, all=T)
me[is.na(me)] <- 0
with(me, summary(areaO - areaE))
with(me, plot(areaO, areaE, pch="."))



o0 <- with(x, aggregate(list(areaO=c(area*antc, area*(1-antc))), data.frame(volume=c(v1ck,v1nck), diameter=c(d1ck, d1nck), siteIndex=sik, species=rep(c("coniferous", "nonConiferous"), each=length(area)), slope=slpk, transportDistance=distk, owner=ownerk), FUN=sum))
e0 <- cbind(y[,grep("^step[0-9]+\\.", names(y), invert=TRUE)], areaE=rowSums(y[,grep("step1\\.", colnames(y))]))
me <- merge(o0, e0, all=T)
me[is.na(me)] <- 0
with(me, summary(areaO - areaE))
with(me, plot(areaO, areaE, pch="."))

o0 <- with(x, aggregate(list(areaO=c(area*antc, area*(1-antc))), data.frame(volume=c(v1ck,v1nck), diameter=c(d1ck, d1nck), siteIndex=sik, species=rep(c("coniferous", "nonConiferous"), each=length(area))), FUN=sum))
e0 <- with(y, aggregate(data.frame(areaE=rowSums(y[,grep("step1\\.", colnames(y))])), data.frame(volume, diameter, siteIndex, species), FUN=sum))
me <- merge(o0, e0, all=T)
me[is.na(me)] <- 0
with(me, summary(areaO - areaE))
with(me, plot(areaO, areaE, pch="."))
