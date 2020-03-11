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

#Use Stem as stemnumber
#dat$stem0 <- pmin(nClasses["stem"], pmax(1, round(nClasses["stem"] - dat$vol0 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))
#dat$stem1 <- pmin(nClasses["stem"], pmax(1, round(nClasses["stem"] - dat$vol1 * nClasses["stem"] / nClasses["vol"] + runif(NROW(dat), -1, 1))))

#Allow move only in one directions
#tt <- dat$vol1 < dat$vol0
#dat$vol1[tt] <- dat$vol0[tt]
#tt <- dat$stem1 > dat$stem0
#dat$stem1[tt] <- dat$stem0[tt]

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
write.table(dat[c("si", "vol0", "vol1", "stem0", "stem1")], "./dat/dNoman.txt", row.names = FALSE)
#
#estimate transition probabilities
cat("data ./dat/dNoman.txt
prior uninformative
si=1
", file="./dat/estiminput.txt")
#
writeLines(sapply(names(nClasses), function(x) {paste(x,paste(1:nClasses[x],collapse = " "))}), "./dat/factors.txt")
#
source("../../EFDMcode/nea/efdmsetuptools.r")
source("../../EFDMcode/nea/efdmestim.r")
pre.estimate("./dat/estiminput.txt", "./dat/factors.txt", "vol stem", "./dat/pNomgmt.RData") #The order need to be "vol stem" and not "stem vol"
#
#Have a look on the Transition Probabilities
tt <- get(load("./dat/pNomgmt.RData"))
rownames(tt) <- apply(expand.grid(sapply(nClasses[c("vol", "stem")], function(x) 1:x)),1,paste,collapse="-") #volClass-stemClass
colnames(tt) <-  rownames(tt)
#colname (volClass-stemClass) gives where you are. The values in the colum
#give you the propability where it will be (given in rowname) in the next
#time step
t1 <- with(dat[dat$vol0==1 & dat$stem0==1,], table(vol1, stem1)); t1/sum(t1)
tt[tt[,1,1]>0, 1, 1]
t1 <- with(dat[dat$vol0==5 & dat$stem0==3,], table(vol1, stem1)); t1/sum(t1)
tt[tt[,rownames(tt) == "3-5",1]>0, rownames(tt) == "3-5", 1]

cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 1
OUTPUTREQUEST_FILENAME 0
", file="./dat/efdminput.txt")


