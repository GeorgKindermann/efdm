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
breaks$owner <- factor(levels(dat$speciesGroup), levels(dat$speciesGroup))

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
  y <- with(x, data.frame(yield, owner, hangneig, rueck, mantyp, seehoehe, speciesGroup0 = speciesGroup, speciesGroup1 = speciesGroup, volume0, volume1, diameter0, diameter1, area = pmin(a0, a1)))
  y <- y[y$area > 0,]
#  tmp <- if(sum(x$a1 > 0) > 0 && nrow(x) > 1) {
  if(nrow(x) > 1) {
    i <- x$aDif < 0
    k <- which(!i)
    rbind(y, do.call(rbind, lapply(which(i), function(j) {
      aggregate(area ~ ., data.frame(x[j, c("yield", "owner", "hangneig", "rueck", "mantyp", "seehoehe")], speciesGroup0 = x$speciesGroup[j], speciesGroup1 = x$speciesGroup[k], volume0 = x$volume0[j], volume1 = x$volume1[k], diameter0 = x$diameter0[j], diameter1 = x$diameter1[k], area = -x$aDif[j]*x$aDif[k]/sum(x$aDif[k]), row.names = NULL), sum)
    })))} else {y}
})
#t2 <- apply(t2$x, 2, function(x) do.call(cbind, x)) #Fuer einen Punkt
t2 <- do.call(rbind, apply(t2$x, 1, function(x) do.call(cbind, x)))
t2 <- aggregate(area ~ ., t2, sum)
env <- efdmTransitionGet(t2, area="area", t0="0", t1="1", breaks=breaks, getArea=TRUE)
env <- efdmTransitionXda(env)
#env <- efdmTransitionSimilar(env, absolute=TRUE, wgt=c(5,10,1,1,1,1,10,7,5), fixed=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE))
#env <- efdmTransitionShifter(env, c(0, 0, 0, 0, 0, 0, 0, 1, 1))
saveRDS(env, file="./dat/transition.RData", compress="xz")


#t1 <- efdmNextArea(state0Area, env, "share")



#Removals during transitions
tt <- with(dat, data.frame(volume0=v0, volume1=v1, diameter0=d0, diameter1=d1, yield=si, speciesGroup, removalTyp))
tt <- cbind(dat[,c("rw","hw","pbfl","ant0","ant1","area","v0","v1")], efdmClassGet(tt, TRUE, breaks))
tt$diameter0[is.na(tt$diameter0)] <- -1 #or addNA
tt$diameter1[is.na(tt$diameter1)] <- -1

t1 <- aggregate(cbind(a0=ant0*area, a1=ant1*area, v0=v0/ifelse(ant0>0, ant0, 1), v1=v1/ifelse(ant1>0, ant1, 1)) ~ rw + hw + pbfl + yield + speciesGroup + removalTyp + volume0 + volume1 + diameter0 + diameter1, tt, sum)

t2 <- aggregate(seq_len(nrow(t1)), list(t1$rw, t1$hw, t1$pbfl), function(z) {
  x <- t1[z,]
  x$aDif <- x$a1 - x$a0
  tmp <- if(sum(x$a1 > 0) > 0) {
    y <- with(x, data.frame(area = pmin(a0, a1), speciesGroup0 = speciesGroup, speciesGroup1 = speciesGroup, volume0, volume1, diameter0, diameter1, removalTyp, v0, v1))
    y <- y[y$area > 0,]
    i <- x$aDif < 0
    rbind(y, do.call(rbind, lapply(which(i), function(j) {
      k <- which(!i)
      data.frame(area = -x$aDif[j]*x$aDif[k]/sum(x$aDif[k]), speciesGroup0 = x$speciesGroup[j], speciesGroup1 = x$speciesGroup[k], volume0 = x$volume0[j], volume1 = x$volume1[k], diameter0 = x$diameter0[j], diameter1 = x$diameter1[k], removalTyp = x$removalTyp[j], v0 = x$v0[j], v1 = x$v1[j])
    })))} else {
          y <- with(x, data.frame(area = a0, speciesGroup0 = speciesGroup, speciesGroup1 = speciesGroup, volume0, volume1=0, diameter0, diameter1=0, removalTyp, v0, v1))
          y[y$area > 0,]
        }
  cbind(yield=x$yield[1], aggregate(cbind(area, v0=v0*area, v1=v1*area) ~ ., tmp, sum), row.names = NULL)
})
t2 <- do.call(rbind, apply(t2$x, 1, function(x) do.call(cbind, x)))
t2 <- aggregate(cbind(area, v0, v1) ~ ., t2, sum)
t2 <- merge(t2, aggregate(cbind(as=area) ~ ., t2[,setdiff(colnames(t2), c("removalTyp","v0","v1"))], sum))
t2$v0 <- t2$v0 / t2$area
t2$v1 <- t2$v1 / t2$area
t2$area <- t2$area / t2$as
t2 <- t2[t2$removalTyp < 4,] #Type 1 - 3 have some removals
t2$harv <- t2$v1 - t2$v0 #As it is single tree based v1 is in all cases 0

i0 <- c("yield", "speciesGroup0", "volume0", "diameter0")
#i1 <- c("yield", "speciesGroup1", "volume1", "diameter1")
#t2 <- data.frame(from=efdmIndexGet(t2[i0], efdmDimGet(i0, breaks)), to=efdmIndexGet(t2[i1], efdmDimGet(i1, breaks)), t2[c("area","removalTyp","harv")])
#t2$removalTyp <- breaks$removalTyp[t2$removalTyp]
#t2 <- reshape(t2, v.names=c("area", "harv"), timevar="removalTyp", idvar=c("from","to"), direction="wide")
#t2[is.na(t2)] <- 0
t2 <- data.frame(from=efdmIndexGet(t2[i0], efdmDimGet(i0, breaks)), t2[c("area","removalTyp","harv")])
t2$removalTyp <- breaks$removalTyp[t2$removalTyp]
t2 <- reshape(t2, v.names=c("area", "harv"), timevar="removalTyp", idvar="from", direction="wide")
t2[is.na(t2)] <- 0

envF <- new.env(parent=emptyenv())
envF$from <- t2$from
#envF$to <- t2$to
envF$share <- as.matrix(t2[,startsWith(names(t2), "area")])
envF$ammount <- as.matrix(t2[,startsWith(names(t2), "harv")])
envF$dimS <- efdmDimGet(i0, breaks)
saveRDS(envF, file="./dat/flow.RData", compress="xz")



state1Area <- efdmNextArea(state0Area, env, "share")


#Flow
t0 <- state0Area
transition <- env
res <- Matrix::sparseVector(0, 1, length(t0))
attr(res, "dimS") <- attr(t0, "dimS")
me <- merge(cbind(efdmIndexInv(t0@i, attr(t0, "dimS")), area=t0@x), cbind(efdmIndexInv(envF$from, envF$dimS), envF$share, envF$ammount))
