source("~/prg/efdm/efdm.r")

dat <- read.table("./dat/owiEfdmBase.txt.bz2", header=T)


dat$removalTyp <- "no"
dat$removalTyp[is.na(dat$d1)] <- "calamity"
dat$removalTyp[!is.na(dat$nutzar) & dat$nutzar=="0"] <- "mortality"
dat$removalTyp[!is.na(dat$nutzar) & dat$nutzar %in% c("1","2","3","4","5","6","7","8")] <- "harvest"
dat$removalTyp <- as.factor(dat$removalTyp)
table(dat$removalTyp)

dat$speciesGroup <- "coniferous"
dat$speciesGroup[dat$bart > 9] <- "nonConiferous"
dat$speciesGroup[dat$bart == 1] <- "spruce"
dat$speciesGroup[dat$bart %in% 4:7] <- "pine"
dat$speciesGroup[dat$bart == 10] <- "beech"
dat$speciesGroup[dat$bart == 11] <- "oak"
dat$speciesGroup <- as.factor(dat$speciesGroup)
table(dat$speciesGroup)


quantileWeighted <- function(x, w, probs = seq(0, 1, 0.25)) {
  x <- aggregate(w[w>0] ~ x[w>0], FUN=sum)
  approxfun(filter(c(0,cumsum(x$w)/sum(x$w)), c(.5,.5), sides=1)[-1], x$x, rule=2)(probs)
}

c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)


#summary(dat$v0)
#summary(dat$v0 / dat$ant0) #Hectare Value
#summary(dat$area * dat$ant0) #Total area


##Define ranges of classes
breaks <- new.env(parent=emptyenv())

#Search for volume classes
with(dat[!is.na(dat$d0) & !is.na(dat$d1),], weighted.mean(v1/ant1 - v0/ant0, (ant0+ant1)*area))
round(with(dat[!is.na(dat$d0) & !is.na(dat$d1),], quantileWeighted(v1/ant1 - v0/ant0, (ant0+ant1)*area, seq(0, 1, 0.05))))
round(with(dat, quantileWeighted(c(v0/ant0, v1/ant1), c(ant0*area, ant1*area), seq(0, 1, 0.05))))
breaks$volume <- c(0.1, seq(25, 1700, 25))

tt <- with(dat, aggregate(cbind(v0*area, ant0*area), list(findInterval(v0/ant0, breaks$volume)), sum))
tt <- data.frame(i=tt[,1], x=tt[,2] / tt[,3])
a <- loess(x ~ i, tt, control = loess.control(surface = "direct"))
tt <- approxfun(tt[,1], tt[,2])(0:length(breaks$volume))
tt[is.na(tt)] <- predict(a, 0:length(breaks$volume))[is.na(tt)]
breaks$stock <- tt
rm(tt, a)

#Search for diameter classes
with(dat[!is.na(dat$d0) & !is.na(dat$d1),], weighted.mean(d1 - d0, (ant0+ant1)*area))
with(dat[!is.na(dat$d0) & !is.na(dat$d1),], quantileWeighted(d1 - d0, (ant0+ant1)*area, seq(0, 1, 0.05)))
with(dat[!is.na(dat$d0) & !is.na(dat$d1),], quantileWeighted(c(d0, d1), c(ant0*area, ant1*area), seq(0, 1, 0.05)))
breaks$diameter <- seq(2, 70, 2)

#Yield classes
(t1 <- with(dat, quantileWeighted(si, (ant0+ant1)*area,  seq(0, 1, 0.05))))
diff(t1)
breaks$yield <- seq(14, 46, 1)

#Species
breaks$speciesGroup <- factor(levels(dat$speciesGroup), levels(dat$speciesGroup))

#Owner type
table(dat$ea)
dat$owner <- factor(dat$ea==1, levels=c(TRUE,FALSE), labels=c("under200ha","other"))
breaks$owner <- factor(levels(dat$owner), levels(dat$owner))

#Slope
table(dat$hangneig)
breaks$hangneig <- c(4, 6, 9)

#Transport distance
summary(dat$rueck)
summary(dat$vorrueck)
breaks$rueck <- c(150, 300, 700)

#Management type
table(dat$ba) #0..Hochwald, 1..HochwaldSchutz, 2..Ausschlagwald
dat$mantyp <- factor(ifelse(dat$ba>2000, 2, ifelse(floor(dat$ba/10) %% 10 != 1, 1, 0)), levels=0:2, labels=c("seedRegen", "seedRegenProtection", "coppice"))
breaks$mantyp <- factor(levels(dat$mantyp), levels(dat$mantyp))

#Altitude
table(dat$seehoehe)
breaks$seehoehe <- c(3, 6, 9, 12, 15, 18)

#removalTyp
breaks$removalTyp <- factor(levels(dat$removalTyp), levels(dat$removalTyp))

saveRDS(breaks, file="./dat/breaks.RData")


#Initial state from where the simulation starts: How much area is in the defined classes
tt <- with(dat[!is.na(dat$d0),], data.frame(volume=v0/ifelse(ant0>0, ant0, 1), diameter=d0, yield=si, speciesGroup, owner, hangneig, rueck, mantyp, seehoehe, area=ant0*area))
state0Area <- efdmStateGet(tt, "area", breaks, TRUE)
saveRDS(state0Area, file="./dat/state0Area.RData", compress="xz")
rm(tt)

#observed t1
tt <- with(dat[!is.na(dat$d1),], data.frame(volume=v1/ifelse(ant1>0, ant1, 1), diameter=d1, yield=si, speciesGroup, owner, hangneig, rueck, mantyp, seehoehe, area=ant1*area))
tt <- efdmStateGet(tt, "area", breaks, TRUE)
saveRDS(tt, file="./dat/state1Area.RData", compress="xz")
rm(tt)


#Transition Probabilities
tt <- with(dat, data.frame(volume0=v0/ifelse(ant0>0, ant0, 1), volume1=v1/ifelse(ant1>0, ant1, 1), diameter0=d0, diameter1=d1, yield=si, speciesGroup, owner, hangneig, rueck, mantyp, seehoehe))
#tt <- cbind(dat[,c("rw","hw","pbfl","ant0","ant1","area","removalTyp")], efdmClassGet(tt, TRUE, breaks))
tt <- cbind(dat[,c("rw","hw","pbfl","ant0","ant1","area","removalTyp")], tt)
tt$diameter0[is.na(tt$diameter0)] <- -1L #or addNA
tt$diameter1[is.na(tt$diameter1)] <- -1L

t1 <- aggregate(cbind(a0=ant0*area, a1=ant1*area) ~ rw + hw + pbfl + yield + owner + hangneig + rueck + mantyp + seehoehe + speciesGroup + removalTyp + volume0 + volume1 + diameter0 + diameter1, tt, sum)

t2 <- aggregate(seq_len(nrow(t1)), list(t1$rw, t1$hw, t1$pbfl), function(z) {
  x <- t1[z,]
  x$aDif <- x$a1 - x$a0
  y <- with(x, data.frame(yield, owner, hangneig, rueck, mantyp, seehoehe, removalTyp, speciesGroup0 = speciesGroup, speciesGroup1 = speciesGroup, volume0, volume1, diameter0, diameter1, area = pmin(a0, a1)))
  y <- y[y$area > 0,]
#  tmp <- if(sum(x$a1 > 0) > 0 && nrow(x) > 1) {
  if(nrow(x) > 1) {
    i <- x$aDif < 0
    k <- which(!i)
    rbind(y, do.call(rbind, lapply(which(i), function(j) {
      aggregate(area ~ ., data.frame(x[j, c("yield", "owner", "hangneig", "rueck", "mantyp", "seehoehe", "removalTyp")], speciesGroup0 = x$speciesGroup[j], speciesGroup1 = x$speciesGroup[k], volume0 = x$volume0[j], volume1 = x$volume1[k], diameter0 = x$diameter0[j], diameter1 = x$diameter1[k], area = -x$aDif[j]*x$aDif[k]/sum(x$aDif[k]), row.names = NULL), sum) })))
  } else {y}
})
#t2 <- apply(t2$x, 2, function(x) do.call(cbind, x)) #Fuer einen Punkt
t2 <- do.call(rbind, apply(t2$x, 1, function(x) do.call(data.frame, x)))
t2 <- aggregate(area ~ ., t2, sum)
env <- efdmTransitionGet(t2[,setdiff(colnames(t2), "removalTyp")], area="area", t0="0", t1="1", breaks=breaks, getArea=TRUE)
#efdmTransitionXda(env)
efdmTransitionSimilar(env, absolute=TRUE, wgt=c(5,99,1,1,1,1,10,7,5), fixed=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE))
#efdmTransitionShifter(env, c(0, 0, 0, 0, 0, 0, 0, 1, 1))
saveRDS(env, file="./dat/transition.RData", compress="xz")

#Removals during transitions
envF <- efdmTransitionFlow(t2[c("area", "volume0", "volume1", "diameter0", "diameter1", "speciesGroup0", "speciesGroup1", "removalTyp", "yield")], breaks=breaks)
efdmTransitionFlowMod(envF)
saveRDS(envF, file="./dat/flow.RData", compress="xz")


