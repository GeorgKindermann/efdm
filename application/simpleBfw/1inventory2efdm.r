set.seed(0)
dir.create("dat")
n <- 1000 #Number of forest stands
dt <- 10 #Timespann between two observations

#Create a simple example dataset
#age0..Age at timestep 0
#age1..Age at timestep 1 (typical age0+dt when there was not clearcut
#area..Area of the plot (hectar)
#drain..Removals and mortality between age0 and age1
#man..Managementtype
#plotId..Plot identification
#si..Siteindex
#vol0..Volume per hectare (area) at timestep 0
#vol1..Volume per hectare (area) at timestep 1
#I added diameter instead of stem-number as age substitution
#But unfortunately you have to name it "stem" (and the other axis "vol")
#stem0..Diamter at timestep 0
#stem1..Diamter at timestep 1
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
} #Truncated normal distribution https://stackoverflow.com/a/19350640/10488504
#Round to given total https://stackoverflow.com/a/55610457/10488504
#x..numeric vector
#target..sum of rounded x, if not given target = round(sum(x), digits)
#digits..number of decimal places
#closest..Make adjustemnt by changing closest number
#ref..reference level to calculate probability of adjustemnt, if ref==NA the probability of an adjustement is equal for all values of x
#random..sholud the adjustment be done stocastic or randomly
roundt <- function(x, target=NA, digits = 0, closest=TRUE, ref=0, random=FALSE, min=-Inf, max=Inf) {
  if(is.na(target)) {target <- round(sum(x), digits)}
  if(all(x == 0)) {
    if(target == 0) {return(x)}
    x <- x + 1
  }
  xr <- round(x, digits)
  bias <- target - sum(xr)
  if(bias == 0) {return(xr)}
  if(is.na(ref)) {wgt <- rep(1, length(x))
  } else {
    if(closest) {wgt <- (x - xr) * sign(target - sum(xr)) + 10^-digits / 2
    } else {wgt <- abs(x-ref)}
  }
  if(bias > 0) {wgt[xr==max] <- 0
  } else {wgt[xr==min] <- 0}
  if(sum(wgt)!=0) {wgt <- wgt / sum(wgt)
  } else {return(xr)}
  if(random) {adj <- table(sample(factor(1:length(x)), size=abs(bias)*10^digits, replace = T, prob=wgt))*sign(bias)*10^-digits
  } else {adj <- diff(c(0,round(cumsum((bias) * wgt), digits)))}
  xr + adj
}
funGwl <- function(age, si) { #Estimate total increment until age for si
  c <- 100*si / (150 * si * (1 - exp((-0.028 + 0.00048*si)*100))^(13.5*si^-0.57))
  c * 150 * si * (1 - exp((-0.028 + 0.00048*si)*age))^(13.5*si^-0.57)
}
funStock <- function(age, si, thin=250) {
  funGwl(age, si) / (1 + age/thin)
}
funUopt <- function(age, si) {funGwl(age, si)/age}
funGetBg <- function(age, si, stock) {
  tt <- funStock(age, si)
  pmin(1, ifelse(tt>0, stock / tt, 1))
}
funIncRed <- function(age, si, bg) {
  bg^pmax(.2, pmin(1, (1 - 2*si/age)))^2
}
funD <- function(age, si) {
  pmax(0, -1 + age * sqrt(si) / 10)
}
funDAdjBg <- function(age, si, bg, d) {
  ifelse(d>0, pmin(2, sqrt(d^2 * funIncRed(age, si, bg) / bg) / d), 1)
}
dat <- data.frame(plotId = 1:n, si=round(runif(n,4,16),1), man=factor("nomgmt", levels = c("nomgmt","thin","clearcut","selectiveHarvest")), area=round(runif(n,0.05,10),2))
t1 <- sapply(dat$si, function(x) optimize(funUopt, interval=c(15,199), maximum=T, si=x)$maximum) #Rotation time to reach GWLmax
t1 <- t1 * rtnorm(n,mean=1,sd=.2,.3,1.7) #Add some nois
dat$age0 <- round(runif(n,0,t1))
dat$age1 <- dat$age0 + dt
dat$vol0 <- funStock(dat$age0, dat$si, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$si/40),.85)
dat$stem0 <- funD(dat$age0, dat$si) * rtnorm(n, 1, .1, .5, 1.5)
dat$stem1 <- dat$stem0 + (funD(dat$age1, dat$si) - funD(dat$age0, dat$si)) * funDAdjBg(dat$age0, dat$si, funGetBg(dat$age0, dat$si, dat$vol0), dat$stem0) * rtnorm(n, 1, .1, .5, 1.5)
t2 <- (funGwl(dat$age1, dat$si) - funGwl(dat$age0, dat$si)) * funIncRed(dat$age0, dat$si, funGetBg(dat$age0, dat$si, dat$vol0)) #increment
t2 <- t2 * rtnorm(n, mean = 1, sd = .1, .5, 1.5) #Add some nois
dat$vol1 <- dat$vol0 + t2
t3 <- floor(runif(n,0,dt))  #Years until managent action
t4 <- (dat$vol0 + (t2*t3/dt)) - funStock(dat$age0+t3, dat$si, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$si/40),0.85) #Thin more in young with high si
dat$man[dat$vol0 + t2 > funStock(dat$age1, dat$si, 1000) * runif(n,0.85,1) & t4 > 0] <- "thin"
dat$man[dat$age1 > t1] <- sample(c("clearcut", "nomgmt"), size=sum(dat$age1 > t1), replace=T, prob=c(.7,.3))
dat$man[t4 > 0 & runif(n) < .3 & funGetBg(dat$age0, dat$si, dat$vol0) > .8] <- "selectiveHarvest"
dat$drain <- 0
tt <- dat$man=="clearcut"
dat$drain[tt] = dat$drain[tt] + dat$vol0[tt] + (t2*t3/dt)[tt]
dat$age1[tt] <- (dt - t3 - floor(runif(n,0,4)))[tt]
dat$vol1[tt] <- (funStock(pmax(0,dat$age1), dat$si, 1000) * runif(n,0.5,1))[tt]
dat$stem1[tt] <- (funD(dat$age1, dat$si) * rtnorm(n, 1, .1, .5, 1.5))[tt]
tt <- dat$man%in%c("thin","selectiveHarvest")
dat$drain[tt] <- dat$drain[tt] + t4[tt]
dat$vol1[tt] <- ((dat$vol0 + (t2*t3/dt)) - t4)[tt]
dat$vol1[tt] <- (dat$vol1 + (funGwl(dat$age1, dat$si) - funGwl(dat$age0+t3, dat$si)) * funIncRed(dat$age0+t3, dat$si, funGetBg(dat$age0+t3, dat$si, dat$vol1)) * rtnorm(n, mean = 1, sd = .1, .5, 1.5))[tt]
dat$stem1[dat$man=="thin"] <- pmax(1, dat$stem1 + rtnorm(n, mean = 2, sd = 1, -2, 4))[dat$man=="thin"]
dat$stem1[dat$man=="selectiveHarvest"] <- pmax(1, dat$stem1 + rtnorm(n, mean = -4, sd = 2, -6, 2))[dat$man=="selectiveHarvest"]
dat$age1[dat$man=="selectiveHarvest"] <- pmax(1, dat$age1 + rtnorm(n, mean = -15, sd = 10, -30, 0))[dat$man=="selectiveHarvest"]
dat$vol0 <- round(dat$vol0, 1)
dat$vol1 <- round(dat$vol1, 1)
dat$drain <- round(dat$drain, 1)
dat$stem0 <- round(dat$stem0, 1)
dat$stem1 <- round(dat$stem1, 1)

#Define ranges of classes
#large volume span to avoid hiting the “roof” of the volume dimension
#probability of growing out of a state for each +- same
#number of volume classes should be equal to the number of time periods (age classes) that it takes to reach its maximum
classVol <- round(funGwl(seq(0, max(c(dat$age0, dat$age1)+dt), dt), weighted.mean(dat$si, dat$area)))
classVol <- classVol[classVol < max(c(dat$vol0, dat$vol1))]
classVol <- round(c(classVol, seq(max(classVol), max(c(dat$vol0, dat$vol1))+2, length.out=2)[-1]))
#classAge <- seq(0,max(c(dat$age0,dat$age1))+dt,dt) #Age classes (should coincides with dt), gk:and max should be considered
classSi <- floor(seq(min(dat$si), max(dat$si)+2, 2)) #Siteindex classes
classStem <- seq(0,max(c(dat$stem0,dat$stem1))+5,5)
nClasses <- c(length(classVol), length(classStem), length(classSi))-1
names(nClasses) <- c("vol", "stem", "si")

#Classify the dataset
dat$cstem0 <- factor(findInterval(dat$stem0, classStem), levels=seq_along(classStem[-1]), ordered=T)
dat$cstem1 <- factor(findInterval(pmax(0,dat$stem1), classStem), levels=seq_along(classStem[-1]), ordered=T)
dat$cs <- factor(findInterval(dat$si, classSi), levels=seq_along(classSi[-1]), ordered=T)
dat$cv0 <- factor(findInterval(dat$vol0, classVol), levels=seq_along(classVol[-1]), ordered=T)
dat$cv1 <- factor(findInterval(dat$vol1, classVol), levels=seq_along(classVol[-1]), ordered=T)

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)

#Levels which have been used
cat(paste("vol", paste(levels(dat$cv0), collapse = ' '))
, paste("stem", paste(levels(dat$cstem0), collapse = ' '))
, paste("si", paste(levels(dat$cs), collapse = ' ')), sep="\n", file="./dat/factors.txt")

#Initial state from where the simulatin starts: How much area is in the defined classes
#I take here timestep 0, but also timestep 1 could be used alternatively
t1 <- with(dat, aggregate(list(area=area), list(stem=cstem0, si=cs, vol=cv0), FUN=sum))
write.table(t1, "./dat/initstate.txt", row.names = FALSE)
rm(t1)

#Activities: transition probabilities depending on different activities which
#are defined here
cat("nomgmt read ../dat/pNomgmt.RData vol stem
clearcut read ../dat/pClearcut.RData vol stem
thin read ../dat/pThin.RData vol stem
selectiveHarvest read ../dat/pSelectiveHarvest.RData vol stem
", file="./dat/activities.txt")

#Propabilities for harvest activities
#proportions of the area in a state-space cell that must be subject to the respective activities within a given simulation period
t1 <- with(dat, aggregate(list(man=area), list(stem=cstem0, si=cs, vol=cv0, man=man), FUN=sum))
t1 <- reshape(t1, timevar = "man", idvar = c("stem", "si", "vol"), direction="wide", sep="")
t1[is.na(t1)] <- 0
colnames(t1) <- gsub("^man.1", "", colnames(t1))
t1[,4:NCOL(t1)] <- t1[,4:NCOL(t1)] / rowSums(t1[,4:NCOL(t1)])
levels(t1[1])
t1 <- merge(expand.grid(lapply(t1[,1:3], levels)), t1, all.x=T)
#Fill up the missing positions of the matrix
t2 <- cbind(apply(t1[,1:3], 2, as.numeric), Reduce(rbind, lapply(colnames(t1)[4:NCOL(t1)], function(x) data.frame(x, wgt=t1[,x]))))
library(MASS)
a <- lda(as.formula(paste0("x ~ ", paste(names(t2)[1:3], collapse=" + "))), data=t2[!is.na(t2$wgt) & t2$wgt>0,], weights=wgt)
t2 <- !complete.cases(t1)
t1[t2,4:NCOL(t1)] <- predict(a, as.data.frame(apply(t1[,1:3], 2, as.numeric)))$posterior[t2,]
#Round
#t1[,4:NCOL(t1)] <- t(apply(t1[,4:NCOL(t1)],  1, function(x) roundt(x=x, target=1, digits=2, closest=F, min=0, max=1)))
write.table(t1, "./dat/actprobs.txt")
#xtabs(nomgmt ~ stem + vol + si, data=t1)
#xtabs(clearcut ~ stem + vol + si, data=t1)
rm(t1)

#Transition Probabilities
source("../../EFDMcode/bfw/transMat.r")
#Clear Cut (Allow also that it is not on position 1,1)
x <- with(dat[dat$man=="clearcut",], aggregate(list(area=area), list(vol1=cv1, stem1=cstem1, vol0=cv0, stem0=cstem0, si=cs), FUN=sum))
prior <- getPriorIJ(nlevels(x[,1]), nlevels(x[,2]), -Inf, -Inf)
weight <- c(1,1)
tt <- transMatWP(x,weight,prior)
save(tt, file="./dat/pClearcut.RData")
#Other Managementforms
for(managmnt in setdiff(levels(dat$man), "clearcut")) {
  x <- with(dat[dat$man==managmnt,], aggregate(list(area=area), list(vol1=cv1, stem1=cstem1, vol0=cv0, stem0=cstem0, si=cs), FUN=sum))
  trans <- getTrans(x, 0.2)
  x <- rbind(x, fillEmptyTarget(x))
  prior <- getPriorLda(x)
  t1 <- transMat2LT(x[0,], prior)
  t1 <- cutLT(t1, trans)
  prior <- getPriorObs(t1)
  weight <- c(1,1) #Use observations by 100% -> reproduce observation
  #weight <- c(0,0) #Use 100% of Prior
  tt <- transMatWP(x,weight,prior)
  substr(managmnt, 1, 1) <- toupper(substr(managmnt, 1, 1))
  save(tt, file=paste0("./dat/p", managmnt, ".RData"))
}


#Output request file
cat("volume by stem
../dat/orVolumeStem.txt
drain by stem si
../dat/orDrainStemSi.txt
area by stem
../dat/orAreaStem.txt
area by stem si
../dat/orAreaStemSi.txt
drain by vol stem si
../dat/orDrainStemSi.txt
volume by vol stem si
../dat/orVolumeVolStemSi.txt
", file="./dat/outputrequests.txt")
#Average Volume per volume and stem class
t1 <- with(dat, aggregate(list(x=c(vol0,vol1)*area, area=area), list(vol=c(cv0,cv1), stem=c(cstem0,cstem1)), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + stem, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","stem")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, stem, volume=round(t,1))), "./dat/orVolumeStem.txt", row.names = FALSE)
#Average Harvest per volume, stem and si class
t1 <- with(dat, aggregate(list(x=drain*area, area=area), list(vol=as.integer(cv0), stem=as.integer(cstem0), si=as.integer(cs)), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + stem + si, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","stem","si")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, stem, si, drain=round(t,1))), "./dat/orDrainStemSi.txt", row.names = FALSE)
#Area per stem
t1 <- data.frame(stem = 1:nClasses["stem"], area=1)
write.table(t1, "./dat/orAreaStem.txt", row.names = FALSE)
#Area per stem and si class
t1 <- expand.grid(sapply(nClasses[c("stem","si")], function(x) 1:x))
t1$area <- 1
write.table(t1, "./dat/orAreaStemSi.txt", row.names = FALSE)
#Average Volume per volume, stem and si class
t1 <- with(dat, aggregate(list(x=c(vol0,vol1)*area, area=area), list(vol=c(cv0,cv1), stem=c(cstem0,cstem1), si=c(cs,cs)), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + stem + si, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","stem","si")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, stem, si, volume=round(t,1))), "./dat/orVolumeVolStemSi.txt", row.names = FALSE)


#Create efdm control file
cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 3
OUTPUTREQUEST_FILENAME ../dat/outputrequests.txt
", file="./dat/efdminput.txt")
