set.seed(0)
dir.create("dat")

#Define ranges of classes
nClasses <- c(10, 7, 1)
names(nClasses) <- c("vol", "age", "si")

#Create a simple example dataset
dat <- expand.grid(sapply(nClasses[c("age", "si")], function(x) 1:x))
dat <- dat[rep(seq_len(nrow(dat)), round(runif(NROW(dat), 3, 7))),]
names(dat) <- c("age0", "si")
dat$area <- 1
dat$age1 <- dat$age0 + round(runif(NROW(dat), 0, 2))
#Make clearcut, breakdown
dat$age1[dat$age1 > nClasses["age"] | runif(NROW(dat)) < .05] <- 1

dat$vol0 <- pmin(nClasses["vol"], pmax(1, round(dat$age0 * nClasses["vol"] / nClasses["age"] + runif(NROW(dat), -1, 1))))
dat$vol1 <- pmin(nClasses["vol"], pmax(1, round(dat$age1 * nClasses["vol"] / nClasses["age"] + runif(NROW(dat), -1, 1))))

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)

#Initial state from where the simulatin starts
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), with(dat, aggregate(list(area=area), list(age=age0, si=si, vol=vol0), FUN=sum)), all.x=T)
#t1 <- t1[!is.na(t1$area),]  #EFDM needs to have also 0 information 
t1$area[is.na(t1$area)] <- 0
write.table(t1, "./dat/initstate.txt", row.names = FALSE)
rm(t1)

#Activities:
cat("nomgmt read ../dat/pNomgmt.RData vol age
", file="./dat/activities.txt")

#Propabilities for harvest activities
t1 <- expand.grid(sapply(nClasses, function(x) 1:x))
t1$nomgmt <- 1
write.table(t1, "./dat/actprobs.txt")
rm(t1)

#Transition Probabilities
write.table(dat[c("si", "vol0", "vol1", "age0", "age1")], "./dat/dNoman.txt", row.names = FALSE)
#Create prior table
#uninformative (everything stays where it is)
#tt <- diag(nrow = (ncol = nClasses["vol"] * nClasses["age"]))
#Take information from the observations
t1 <- aggregate(area ~ vol0 + vol1 + age0 + age1, data=dat, FUN=sum)
t2 <- expand.grid(sapply(nClasses[c("vol","age")], function(x) 1:x))
names(t2) <- paste0(names(t2),0)
me <- merge(t2, t1, all.x=T)
me$vol1[is.na(me$vol1)] <- me$vol0[is.na(me$vol1)] #No info -> stay there
me$age1[is.na(me$age1)] <- me$age0[is.na(me$age1)]
me$area[is.na(me$area)] <- 1
t2$vol1 <- t2$vol0
t2$age1 <- t2$age0
t2$area <- 0
me <- rbind(me, t2)
tt <- prop.table(xtabs(area~I(age1*nClasses["vol"]+vol1)+I(age0*nClasses["vol"]+vol0), me,addNA=T),2)
write.table(tt, "./dat/dPrior.txt", row.names=F, col.names=F)
#
#estimate transition probabilities
cat("data ./dat/dNoman.txt
prior ./dat/dPrior.txt
si=1
", file="./dat/estiminput.txt")
#
writeLines(sapply(names(nClasses), function(x) {paste(x,paste(1:nClasses[x],collapse = " "))}), "./dat/factors.txt")
#
source("../../EFDMcode/v2.0/hackfunctions.r")
source("../../EFDMcode/v2.0/efdmestim.r")
pre.estimate("./dat/estiminput.txt", "./dat/factors.txt", "vol age", "./dat/pNomgmt.RData") #The order need to be "vol age" and not "age vol"
#
#Have a look on the Transition Probabilities
tt <- get(load("./dat/pNomgmt.RData"))
rownames(tt) <- apply(expand.grid(sapply(nClasses[c("vol", "age")], function(x) 1:x)),1,paste,collapse="-") #volClass-ageClass
colnames(tt) <-  rownames(tt)
#colname (volClass-ageClass) gives where you are. The values in the colum
#give you the propability where it will be (given in rowname) in the next
#time step
t1 <- with(dat[dat$vol0==1 & dat$age0==1,], table(vol1, age1)); t1/sum(t1)
tt[tt[,1,1]>0, 1, 1]
t1 <- with(dat[dat$vol0==7 & dat$age0==5,], table(vol1, age1)); t1/sum(t1)
tt[tt[,rownames(tt) == "7-5",1]>0, rownames(tt) == "7-5", 1]

cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 1
OUTPUTREQUEST_FILENAME 0
", file="./dat/efdminput.txt")


