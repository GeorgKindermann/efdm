set.seed(0)
dir.create("dat")

#Define ranges of classes
nClasses <- c(10, 7, 3)
names(nClasses) <- c("vol", "stem", "si")

#Create a simple example dataset
dat <- expand.grid(sapply(nClasses[c("vol", "si")], function(x) 1:x))
dat <- dat[rep(seq_len(nrow(dat)), round(runif(NROW(dat), 1, 7))),]
names(dat) <- c("vol0", "si")
dat$vol1 <- pmin(nClasses["vol"], pmax(1, dat$vol0 + round(runif(NROW(dat), -1, 2+dat$si))))
dat$area <- 1
#Make a breakdown
dat$brk <- FALSE
dat$brk[sample(seq_len(NROW(dat)), round(.1*NROW(dat)))] <- TRUE
dat$vol1[dat$brk] <- pmin(nClasses["vol"], pmax(1, dat$vol1[dat$brk] - ceiling(runif(sum(dat$brk), 1, dat$vol1[dat$brk]))))

#Use Stem as diameter
dat$stem0 <- pmin(nClasses["stem"], pmax(1, round(dat$vol0 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))
dat$stem1 <- pmin(nClasses["stem"], pmax(1, round(dat$vol1 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)

dat$vol0 <- factor(dat$vol0, levels=seq_len(nClasses["vol"]), ordered=TRUE)
dat$vol1 <- factor(dat$vol1, levels=seq_len(nClasses["vol"]), ordered=TRUE)
dat$stem0 <- factor(dat$stem0, levels=seq_len(nClasses["stem"]), ordered=TRUE)
dat$stem1 <- factor(dat$stem1, levels=seq_len(nClasses["stem"]), ordered=TRUE)
dat$si <- factor(dat$si, levels=seq_len(nClasses["si"]), ordered=TRUE)

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

writeLines(sapply(names(nClasses), function(x) {paste(x,paste(1:nClasses[x],collapse = " "))}), "./dat/factors.txt")

#Transition Probabilities: Development without management
source("../../EFDMcode/bfw/transMat.r")
x <- aggregate(area ~ vol1 + stem1 + vol0 + stem0 + si, data=dat, FUN=sum)
#prior <- getPriorObs(x)
prior <- getPriorRF(x)
#prior <- getPriorSvm(x)
#prior <- getPriorIJ(nlevels(dat$vol1),nlevels(dat$stem1),1,1)
weight <- c(1,1) #Use observations by 100% -> reproduce observation
#weight <- c(0,0) #Use 100% of Prior
tt <- transMatWP(x,weight,prior)
save(tt, file="./dat/pNomgmt.RData")

cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 1
OUTPUTREQUEST_FILENAME 0
", file="./dat/efdminput.txt")


