source("efdm.r")

#Read Inventory data
dat <- read.csv("./dat/inventoryData.csv.bz2")

dt <- max(dat$age1 - dat$age0)
n <- nrow(dat)

#Define ranges of classes
#The size of a growth classes (volumen, age, diameter) should be maximum the
#change in one growht period, so that all leave the class in one period
#but it shuld not be too small to have gaps on the transition path
breaks <- new.env(parent=emptyenv())
breaks$age <- seq(10,max(c(dat$age0,dat$age1)),dt)
breaks$volume <- seq(30,max(c(dat$volume0,dat$volume1)),30)
breaks$diameter <- seq(3,max(c(dat$diameter0,dat$diameter1)),3)
breaks$yield <- seq(ceiling(min(dat$yieldC,dat$yieldNC)),max(dat$yieldC,dat$yieldNC),1)
breaks$temperature <- seq(ceiling(min(dat$temperature)),max(dat$temperature),1)
breaks$species <- 1
saveRDS(breaks, file="./dat/breaks.RData")

#Initial state from where the simulation starts: How much area is in the defined classes
#I take here timestep 0, but also timestep 1 could be used alternatively
#Maximum size .Machine$integer.max bzw. 2^53
tt <- efdmClassGet(dat[,c("volume0","diameter0","temperature","yieldC","yieldNC")], TRUE, breaks)
tt <- cbind(tt[,1:3], yield=c(tt$yieldC, tt$yieldNC), species=rep(0:1, each=nrow(tt)), area=c(dat[,"area"] * dat$shareC0, dat[,"area"] * (1-dat$shareC0)))
t0 <- efdmDimGet(colnames(tt)[-ncol(tt)],  breaks)
if(TRUE) {   #Sparse array
  tt <- aggregate(tt["area"], list(i=efdmIndexGet(tt[-ncol(tt)], t0)), sum)
  state0Area <- Matrix::sparseVector(tt$area, tt$i, prod(t0))
} else {     #Dense array
  state0Area <- as.vector(xtabs(area ~ ., data.frame(tt["area"], Map(factor, tt[-ncol(tt)], lapply(t0-1, seq, from=0)))))
}
#object.size(state0Area)
attr(state0Area, "dimS") <- t0
saveRDS(state0Area, file="./dat/state0Area.RData", compress="xz")
rm(tt, t0)

#Transition Probabilities
transition <- new.env(parent=emptyenv())
#
i <- c("area","stock","harvest","residuals")
tmp <- aggregate(cbind(areaC0=shareC0*area, areaC1=shareC1*area, area, stock=(dat$volume0+dat$volume1)/2*area, harvest=harvest*area, residuals=residuals*area) ~ volume0 + volume1 + diameter0 + diameter1 + yieldC + yieldNC, dat, sum)
t0 <- pmin(1, tmp$areaC1/tmp$areaC0)
t0[is.na(t0)] <- 0
t1 <- pmin(1, (tmp$area-tmp$areaC1)/(tmp$area-tmp$areaC0))
t1[is.na(t1)] <- 0
j <- grep("^(diameter|volume)", colnames(tmp))
k <- tmp[match(i[-1], colnames(tmp))] / tmp$area
tmp  <- rbind(data.frame(tmp[j], yield0=tmp$yieldC, yield1=tmp$yieldC, area = tmp$areaC0 * t0, species0 = 0, species1 = 0, k)
, data.frame(tmp[j], yield0=tmp$yieldC, yield1=tmp$yieldNC, area = tmp$areaC0 * (1-t0), species0 = 0, species1 = 1, k)
, data.frame(tmp[j], yield0=tmp$yieldNC, yield1=tmp$yieldNC, area = (tmp$area - tmp$areaC0) * t1, species0 = 1, species1 = 1, k)
, data.frame(tmp[j], yield0=tmp$yieldNC, yield1=tmp$yieldC, area = (tmp$area - tmp$areaC0) * (1-t1), species0 = 1, species1 = 0, k))
tmp <- tmp[tmp$area > 0,]
tmp[i[-1]] <- tmp[i[-1]] * tmp[,i[1]]
rm(j,k,t0,t1)
#
t0 <- efdmDimGet(grep("[01]$", colnames(tmp), value=TRUE),  breaks)
transition$dimS <- t0
tmp <- aggregate(. ~ from + to, data.frame(from = efdmIndexGet(efdmClassGet(tmp[,endsWith(colnames(tmp), "0")], TRUE, breaks), t0)
 , to = efdmIndexGet(efdmClassGet(tmp[,endsWith(colnames(tmp), "1")], TRUE, breaks), t0)
 , tmp[match(i, colnames(tmp))]), sum)
#
tt <- aggregate(cbind(area,stock) ~ idx, cbind(tmp[c("area","stock")], idx=c(tmp$from, tmp$to)), sum)
tt$stock <- tt$stock / tt$area
if(TRUE) { #Sparse
  state01Stock <- Matrix::sparseVector(tt$stock, tt$idx, prod(t0))
} else {   #Dense
  state01Stock <- "[<-"(numeric(prod(t0)), tt$idx, tt$stock)
}
#object.size(state01Stock)
attr(state01Stock, "dimS") <- t0
saveRDS(state01Stock, file="./dat/state01Stock.RData", compress="xz")
#
tmp[i[-1]] <- tmp[i[-1]] / tmp[,i[1]]
tt <- aggregate(area~from, tmp, sum)
tmp$share <- tmp$area / tt$area[findInterval(tmp$from, tt$from)]
tmp <- tmp[,c("from","to","share","harvest","residuals")]
tmp <- tmp[order(tmp$to),]
rm(tt, i)
#
i <- tmp$share == 1 #All to one
transition$allToOne <- tmp[i,c("from","to","harvest","residuals")]
tmp <- tmp[!i,]
i <- aggregate(cbind(n=from) ~ to, tmp, length)
i <- i$n[match(tmp$to, i$to)] < 5 #N split in some and many
transition$fromSome <- split(tmp[i,], ave(tmp$to[i], tmp$to[i], FUN=seq_along))
transition$fromMany <- do.call(rbind, unname(lapply(split(tmp[!i,], tmp$to[!i]), function(x) data.frame(I(list(x$from)), x$to[1], unname(rbind(as.list(x[,-1:-2])))))))
colnames(transition$fromMany) <- colnames(tmp)
attributes(transition$fromMany$from) <- NULL
#Maybe add also possibility to store it in a dense matrix
#object.size(tmp)
#object.size(transition$fromSome)
#object.size(transition$fromMany)
saveRDS(transition, file="./dat/transition.RData", compress="xz")
rm(tmp, i, t0)
