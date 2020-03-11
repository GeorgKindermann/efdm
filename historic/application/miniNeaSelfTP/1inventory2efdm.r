set.seed(0)
dir.create("dat")

#Define ranges of classes
nClasses <- c(10, 7, 1)
names(nClasses) <- c("vol", "stem", "si")

#Create a simple example dataset
dat <- expand.grid(sapply(nClasses[c("vol", "si")], function(x) 1:x))
dat <- dat[rep(seq_len(nrow(dat)), round(runif(NROW(dat), 1, 7))),]
names(dat) <- c("vol0", "si")
dat$vol1 <- pmin(nClasses["vol"], pmax(1, dat$vol0 + round(runif(NROW(dat), -1, 2))))
dat$area <- 1
#Make a breakdown
dat$vol1[sample(seq_len(NROW(dat)), round(.1*NROW(dat)))] <- 1

#Use Stem as diameter
dat$stem0 <- pmin(nClasses["stem"], pmax(1, round(dat$vol0 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))
dat$stem1 <- pmin(nClasses["stem"], pmax(1, round(dat$vol1 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)

#Initial state from where the simulatin starts
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), with(dat, aggregate(list(area=area), list(stem=stem0, si=si, vol=vol0), FUN=sum)), all.x=T)
t1$area[is.na(t1$area)] <- 0
write.table(t1, "./dat/initstate.txt", row.names = FALSE)
rm(t1)

#Activities:
cat("nomgmt read ../dat/pNomgmt.RData vol stem
", file="./dat/activities.txt")

#Propabilities for harvest activities
t1 <- expand.grid(sapply(nClasses, function(x) 1:x))
t1$nomgmt <- 1
write.table(t1, "./dat/actprobs.txt")
rm(t1)

#Transition Probabilities: Development without management
#Create it by my own
t1 <- aggregate(area ~ vol0 + vol1 + stem0 + stem1 + si, data=dat, FUN=sum)
t2 <- expand.grid(sapply(nClasses, function(x) 1:x))
names(t2)[1:2] <- paste0(names(t2)[1:2],0)
me <- merge(t2, t1, all.x=T)
me$vol1[is.na(me$vol1)] <- me$vol0[is.na(me$vol1)]
me$stem1[is.na(me$stem1)] <- me$stem0[is.na(me$stem1)]
me$area[is.na(me$area)] <- 1
t2$vol1 <- t2$vol0
t2$stem1 <- t2$stem0
t2$area <- 0
me <- rbind(me, t2)
t1 <- array(0, dim = c(nClasses["vol"]*nClasses["stem"],nClasses["vol"]*nClasses["stem"],nClasses["si"]))
for(i in seq_len(nClasses["si"])) {
  t1[,,i] <- prop.table(xtabs(area~I(stem1*nClasses["vol"]+vol1)+I(stem0*nClasses["vol"]+vol0), me[me$si==i,],addNA=T),2)
}
save(t1, file="./dat/pNomgmt.RData")

#Levels of the factors
writeLines(sapply(names(nClasses), function(x) {paste(x,paste(1:nClasses[x],collapse = " "))}), "./dat/factors.txt")

cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 1
OUTPUTREQUEST_FILENAME 0
", file="./dat/efdminput.txt")


