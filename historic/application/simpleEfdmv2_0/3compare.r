x <- read.csv("./dat/inventoryData.csv.bz2")

##Compare area
y <- read.table("./res/rawoutput.txt")

#Timestep 0 by management types
t1 <- with(x, aggregate(list(oArea=area), list(vol=cv0, age=ca0, si=cs, man=man), FUN=sum))
t1 <- reshape(t1, v.names = "oArea", timevar = "man", idvar = c("vol", "age", "si"), direction = "wide")
t2 <- with(y, aggregate(list(eArea.clearcut=step0.clearcut, eArea.nomgmt=step0.nomgmt, eArea.thin=step0.thin), list(vol=vol, age=age, si=si), FUN=sum))

me <- merge(t1, t2, all=T)
me[is.na(me)] <- 0

with(me, plot(oArea.clearcut, eArea.clearcut))
abline(0,1)
with(me, plot(oArea.nomgmt, eArea.nomgmt))
abline(0,1)
with(me, plot(oArea.thin, eArea.thin))
abline(0,1)


#Timestep 1
t1 <- with(x, aggregate(list(oArea=area), list(vol=cv1, age=ca1, si=cs), FUN=sum))
t2 <- with(y, aggregate(list(eArea=step1.clearcut+step1.nomgmt+step1.thin), list(vol=vol, age=age, si=si), FUN=sum))

me <- merge(t1, t2, all=T)
me[is.na(me)] <- 0

with(me, plot(oArea, eArea))
abline(0,1)

with(me, boxplot(oArea - eArea ~ vol))
with(me, boxplot(oArea - eArea ~ age))
with(me, boxplot(oArea - eArea ~ si))

for(si in unique(me$si)) {
  tt <- matrix(0, nrow = max(me$age), ncol = max(me$vol))
  i <- me$si==si
  tt[cbind(me$age[i], me$vol[i])] <- me$oArea[i] - me$eArea[i]
  image(1:max(me$age), 1:max(me$vol), tt, xlab="Age Class", ylab="Volume Class", main=paste0("Observed - Estimated Area, SI=",si))
  with(me[me$si==si,], text(age, vol, round(oArea - eArea)))
  readline()
}


##Compare harvest
y <- read.table("./res/drain_by_vol_age_si.txt")

t1 <- with(x, aggregate(list(oDrain=drain*area), list(vol=cv0, age=ca0, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oDrain, step0))
abline(0,1)


##Compare Stock
y <- read.table("./res/volume_by_vol_age_si.txt")

#Timestep 0
sum(x$vol0*x$area)
sum(y$step0)

t1 <- with(x, aggregate(list(oVol=vol0*area), list(vol=cv0, age=ca0, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oVol, step0))
abline(0,1)

#Timestep 1
sum(x$vol1*x$area)
sum(y$step1)

t1 <- with(x, aggregate(list(oVol=vol1*area), list(vol=cv1, age=ca1, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oVol, step1))
abline(0,1)

with(me, plot(step1, oVol - step1))

with(me, boxplot(oVol - step1 ~ vol))
with(me, boxplot(oVol - step1 ~ age))
with(me, boxplot(oVol - step1 ~ si))

for(si in unique(me$si)) {
  tt <- matrix(0, nrow = max(me$age), ncol = max(me$vol))
  i <- me$si==si
  tt[cbind(me$age[i], me$vol[i])] <- me$oVol[i] - me$step1[i]
  image(1:max(me$age), 1:max(me$vol), tt, xlab="Age Class", ylab="Volume Class", main=paste0("Observed - Estimated Stock, SI=",si))
  with(me[me$si==si,], text(age, vol, round(oVol - step1), cex=0.7))
  readline()
}
