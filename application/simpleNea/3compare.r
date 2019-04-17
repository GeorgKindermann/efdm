x <- read.csv("./dat/inventoryData.csv.bz2")

#Recode management activities
tt <- !x$man %in% c("nomgmt", "clearcut")
x$man <- as.character(x$man)
x$man[tt] <- paste0("man", x$dMan[tt])
x$man <- as.factor(x$man)

##Compare area
y <- read.table("./res/rawoutput.txt")

#Timestep 0 by management types
t1 <- with(x, aggregate(list(oArea=area), list(vol=cv0, stem=cstem0, si=cs, man=man), FUN=sum))
t1 <- reshape(t1, v.names = "oArea", timevar = "man", idvar = c("vol", "stem", "si"), direction = "wide")
t2 <- aggregate(y[grep("step0", names(y), value=T)], list(vol=y$vol, stem=y$stem, si=y$si), FUN=sum)

me <- merge(t1, t2, all=T)
me[is.na(me)] <- 0

sum(me[grep("oArea", names(me), value=T)])
sum(me[grep("step0", names(me), value=T)])

for(i in sub("oArea.", "", grep("oArea", names(me), value=T))) {
  plot(me[,paste0("oArea.",i)], me[,paste0("step0.",i)], main=i, xlab="Observed Area", ylab="Model Area")
  abline(0,1)
  readline()
}

#Timestep 1
t1 <- with(x, aggregate(list(oArea=area), list(vol=cv1, stem=cstem1, si=cs), FUN=sum))
t2 <- aggregate(list(eArea=rowSums(y[grep("step1", names(y), value=T)])), list(vol=y$vol, stem=y$stem, si=y$si), FUN=sum)

me <- merge(t1, t2, all=T)
me[is.na(me)] <- 0

sum(me$oArea)
sum(me$eArea)

with(me, plot(oArea, eArea))
abline(0,1)

with(me, boxplot(oArea - eArea ~ vol))
with(me, boxplot(oArea - eArea ~ stem))
with(me, boxplot(oArea - eArea ~ si))

for(si in unique(me$si)) {
  tt <- matrix(0, nrow = max(me$stem), ncol = max(me$vol))
  i <- me$si==si
  tt[cbind(me$stem[i], me$vol[i])] <- me$oArea[i] - me$eArea[i]
  image(1:max(me$stem), 1:max(me$vol), tt, xlab="Stem Class", ylab="Volume Class", main=paste0("Observed - Estimated Area, SI=",si))
  with(me[me$si==si,], text(stem, vol, round(oArea - eArea)))
  readline()
}


##Compare harvest
y <- read.table("./res/drain_by_vol_stem_si.txt")

t1 <- with(x, aggregate(list(oDrain=drain*area), list(vol=cv0, stem=cstem0, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oDrain, step0))
abline(0,1)


##Compare Stock
y <- read.table("./res/volume_by_vol_stem_si.txt")

#Timestep 0
sum(x$vol0*x$area)
sum(y$step0)

t1 <- with(x, aggregate(list(oVol=vol0*area), list(vol=cv0, stem=cstem0, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oVol, step0))
abline(0,1)

#Timestep 1
sum(x$vol1*x$area)
sum(y$step1)

t1 <- with(x, aggregate(list(oVol=vol1*area), list(vol=cv1, stem=cstem1, si=cs), FUN=sum))
me <- merge(t1, y, all=T)
me[is.na(me)] <- 0

with(me, plot(oVol, step1))
abline(0,1)

with(me, plot(step1, oVol - step1))

with(me, boxplot(oVol - step1 ~ vol))
with(me, boxplot(oVol - step1 ~ stem))
with(me, boxplot(oVol - step1 ~ si))

for(si in unique(me$si)) {
  tt <- matrix(0, nrow = max(me$stem), ncol = max(me$vol))
  i <- me$si==si
  tt[cbind(me$stem[i], me$vol[i])] <- me$oVol[i] - me$step1[i]
  image(1:max(me$stem), 1:max(me$vol), tt, xlab="Stem Class", ylab="Volume Class", main=paste0("Observed - Estimated Stock, SI=",si))
  with(me[me$si==si,], text(stem, vol, round(oVol - step1), cex=0.7))
  readline()
}
