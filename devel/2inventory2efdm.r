#Read Inventory data
dat <- read.csv("./dat/inventoryData.csv.bz2")

dt <- max(dat$age1 - dat$age0)
n <- nrow(dat)

c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

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
saveRDS(breaks, file="./dat/breaks.RData")

#Initial state from where the simulation starts: How much area is in the defined classes
#I take here timestep 0, but also timestep 1 could be used alternatively
#Maximum size .Machine$integer.max bzw. 2^53
t0 <- tt <- c("volume","diameter","temperature")
t0[1:2]  <- paste0(t0[1:2], 0)
key <- as.data.frame(lapply(setNames(seq_along(tt),tt), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- rep("yield", 2)
t0 <- paste0(tt, c("C", "NC"))
keyY <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- cbind(key, yield=c(keyY$yieldC, keyY$yieldNC), species=factor(rep(0:1, each=n), levels=0:1), area=c(dat[,"area"] * dat$shareC0, dat[,"area"] * (1-dat$shareC0)))
t0 <- cumprod(c(1,sapply(tt[,-ncol(tt)], nlevels)))
names(t0) <- names(t0)[-1]
#
#Sparse array
#library(slam)
#library(spray)
tt$idx <- as.vector(1 + matrix(unlist(lapply(tt[,-ncol(tt)], function(f) as.numeric(levels(f))[f])), ncol=length(t0)-1) %*% t0[-length(t0)])
tt <- aggregate(area ~ idx, tt, sum)
state0Area <- Matrix::sparseVector(tt$area, tt$idx, unname(t0[length(t0)]))
attr(state0Area, "idxMul") <- t0[-length(t0)]
saveRDS(list(idxMul=attr(state0Area, "idxMul"), state0Area@x, state0Area@i, state0Area@length), file="./dat/state0Area.RData", compress="xz")
#or Dense array
#state0Area <- as.vector(xtabs(area ~ ., tt))
#attr(state0Area, "idxMul") <- t0[-length(t0)]
#saveRDS(state0Area, file="./dat/state0Area.RData", compress="xz")
#object.size(stateInit)
rm(tt)

#Transition Probabilities
transition <- new.env(parent=emptyenv())
t0 <- tt <- c(rep(c("volume","diameter"),2))
t0[1:2]  <- paste0(t0[1:2], 0)
t0[3:4]  <- paste0(t0[3:4], 1)
key <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- rep("yield", 2)
t0 <- paste0(tt, c("C", "NC"))
keyY <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tmp <- aggregate(cbind(areaC0=shareC0*area, areaC1=shareC1*area, area, stock=(dat$volume0+dat$volume1)/2*area, harvest=harvest*area, residuals=residuals*area) ~ . + volume0 + volume1, cbind(key, keyY, dat[c("shareC0","shareC1","area","harvest","residuals")]), sum)
tmp <- as.data.frame(lapply(tmp, function(x) {if(is.factor(x)) as.numeric(levels(x))[x] else x}))
t0 <- pmin(1, tmp$areaC1/tmp$areaC0)
t0[is.na(t0)] <- 0
t1 <- pmin(1, (tmp$area-tmp$areaC1)/(tmp$area-tmp$areaC0))
t1[is.na(t1)] <- 0
i <- grep("^(diameter|volume)", colnames(tmp))
j <- tmp[match(c("stock","harvest","residuals"), colnames(tmp))]/tmp$area
tmp  <- rbind(data.frame(tmp[i], species0 = 0, species1 = 0, yield0=tmp$yieldC, yield1=tmp$yieldC, area = tmp$areaC0 * t0, j)
, data.frame(tmp[i], species0 = 0, species1 = 1, yield0=tmp$yieldC, yield1=tmp$yieldNC, area = tmp$areaC0 * (1-t0), j)
, data.frame(tmp[i], species0 = 1, species1 = 1, yield0=tmp$yieldNC, yield1=tmp$yieldNC, area = (tmp$area - tmp$areaC0) * t1, j)
, data.frame(tmp[i], species0 = 1, species1 = 0, yield0=tmp$yieldNC, yield1=tmp$yieldC, area = (tmp$area - tmp$areaC0) * (1-t1), j))
rm(j)
i <- match(c("stock","harvest","residuals"), colnames(tmp))
tmp[i] <- tmp[i] * tmp$area
tmp <- aggregate(cbind(area, stock, harvest, residuals) ~ .,tmp[tmp$area>0,],sum)
tt <- lapply(paste0("^", names(attr(state0Area, "idxMul"))), grep, colnames(tmp))
i <- which(lapply(tt, length)==0)
transition$static <- rowSums(expand.grid(lapply(i, function(i) seq(0, by=attr(state0Area, "idxMul")[i], length.out=(c(attr(state0Area, "idxMul")[-1], length(state0Area)) / attr(state0Area, "idxMul"))[i]))))
for(i in seq_len(length(tt))) {
  tmp[,tt[[i]]] <- tmp[,tt[[i]]] * attr(state0Area, "idxMul")[i]
}
i <- match(c("area","stock","harvest","residuals"), colnames(tmp))
tmp <- data.frame(from = rowSums(tmp[,grep("0$", colnames(tmp))])+1
, to = rowSums(tmp[,grep("1$", colnames(tmp))])+1
, tmp[i])
tt <- aggregate(cbind(area,stock) ~ idx, cbind(tmp[c("area","stock")], idx=c(tmp$from, tmp$to)), sum)
tt$stock <- tt$stock / tt$area
#Sparse
state01Stock <- Matrix::sparseVector(tt$stock, tt$idx, length(state0Area))
attr(state01Stock, "idxMul") <- attr(state0Area, "idxMul")
saveRDS(list(idxMul=attr(state01Stock, "idxMul"), state01Stock@x, state01Stock@i, state01Stock@length), file="./dat/state01Stock.RData", compress="xz")
#or dense
#state01Stock <- "[<-"(numeric(length(state0Area)), tt$idx, tt$stock)
#attr(state01Stock, "idxMul") <- attr(state0Area, "idxMul")
#saveRDS(state01Stock, file="./dat/state01Stock.RData", compress="xz")
tmp$stock <- NULL
tt <- aggregate(area~from, tmp, sum)
tmp$share <- tmp$area / tt$area[findInterval(tmp$from, tt$from)]
i <- match(c("harvest","residuals"), colnames(tmp))
tmp[i] <- tmp[i] / tmp$area
tmp$area <- NULL
tmp <- tmp[,c("from","to","share","harvest","residuals")]
tmp <- tmp[order(tmp$to),]
#
i <- tmp$share == 1 #All to one
transition$allToOne <- tmp[i,c("from","to","harvest","residuals")]
tmp <- tmp[!i,]
i <- aggregate(cbind(n=from) ~ to, tmp, length)
i <- i$n[match(tmp$to, i$to)] < 5 #N split in some and many
transition$fromSome <- split(tmp[i,], ave(tmp$to[i], tmp$to[i], FUN=seq_along))
transition$fromMany <- do.call(rbind, unname(lapply(split(tmp[!i,], tmp$to[!i]), function(x) data.frame(to=x$to[1], fs=I(list(unname(as.list(x[-2]))))))))
#Maybe add also possibility to store it in a dense matrix
#object.size(tmp)
#object.size(transition$toSome)
#object.size(transition$toMany)
rm(tmp)
saveRDS(transition, file="./dat/transition.RData", compress="xz")
