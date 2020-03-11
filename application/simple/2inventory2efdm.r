source("../../efdm.r")

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
tt <- cbind(dat[,c("volume0","diameter0","temperature")], yield=c(dat$yieldC, dat$yieldNC), species=rep(0:1, each=nrow(dat)), area=c(dat[,"area"] * dat$shareC0, dat[,"area"] * (1-dat$shareC0)))
colnames(tt) <- sub("[01]$", "", colnames(tt))
state0Area <- efdmStateGet(tt, "area", breaks, TRUE)
#object.size(state0Area)
saveRDS(state0Area, file="./dat/state0Area.RData", compress="xz")
rm(tt)

#Transition Probabilities
i <- c("area","harvest","residuals")
tmp <- aggregate(cbind(areaC0=shareC0*area, areaC1=shareC1*area, area, harvest=harvest*area, residuals=residuals*area) ~ volume0 + volume1 + diameter0 + diameter1 + yieldC + yieldNC, dat, sum)
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
transition <- efdmTransitionGet(tmp, area="area", flow=c("harvest","residuals"), t0="0", t1="1", breaks=breaks)
saveRDS(transition, file="./dat/transition.RData", compress="xz")
rm(i,j,k,t0,t1)
#
state01Stock <- c(0, breaks$volume) + c(diff(c(0, breaks$volume))/2, 0)
tt <- with(tmp, data.frame(area, volume = c(volume0, volume1)))
t0 <- efdmStateGet(tt, "area", breaks, FALSE)
tt$area <- tt$area * tt$volume
tt <- efdmStateGet(tt, "area", breaks, FALSE) / t0
i <- is.finite(tt)
state01Stock[i] <- tt[i]
attr(state01Stock, "dimS") <- attr(t0, "dimS")
saveRDS(state01Stock, file="./dat/state01Stock.RData", compress="xz")
rm(tt, t0, i, tmp)
