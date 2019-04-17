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
dat$cstem0 <- findInterval(dat$stem0, classStem)
dat$cstem1 <- findInterval(pmax(0,dat$stem1), classStem)
dat$cs <- findInterval(dat$si, classSi)
dat$cv0 <- findInterval(dat$vol0, classVol)
dat$cv1 <- findInterval(dat$vol1, classVol)

#Find distinct Management forms
tt <- with(dat[dat$man %in% c("thin", "selectiveHarvest"),], aggregate(list(area=area), list(cstem0=cstem0, cv0=cv0, cstem1=cstem1, cv1=cv1), FUN = sum))
tt <- tt[with(tt, order(cstem0, cv0, -area)),]
tt$dMan <- with(tt, ave(cstem0, cstem0, cv0, FUN=seq_along))
dat <- merge(dat, tt[c("cstem0", "cv0", "cstem1", "cv1", "dMan")], all.x=T)
dat$dMan[!dat$man %in% c("thin", "selectiveHarvest")] <- 0

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)


#Initial state from where the simulatin starts: How much area is in the defined classes
#I take here timestep 0, but also timestep 1 could be used alternatively
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), with(dat, aggregate(list(area=area), list(stem=cstem0, si=cs, vol=cv0), FUN=sum)), all.x=T)
t1$area[is.na(t1$area)] <- 0
write.table(t1, "./dat/initstate.txt", row.names = FALSE)
rm(t1)

#Activities: transition probabilities depending on different activities which
#are defined here
t1 <- "\n"
for(i in unique(tt$dMan)) {
  t1 <- paste0(t1, "man", i, " read ../dat/pMan", i, ".RData vol stem\n")}
cat("nomgmt read ../dat/pNomgmt.RData vol stem
clearcut read ../dat/pClearcut.txt vol stem", t1, file="./dat/activities.txt")
#no management could be estimated on the fly with "nomgmt estimate estiminput.txt vol age"

#Propabilities for harvest activities
#proportions of the area in a state-space cell that must be subject to the respective activities within a given simulation period
t1 <- with(dat, aggregate(list(all=area, nomgmt=area*(man=="nomgmt"), clearcut=area*(man=="clearcut")), list(stem=cstem0, si=cs, vol=cv0), FUN=sum))
t2 <- with(dat[dat$dMan>0,], aggregate(list(man=area), list(stem=cstem0, si=cs, vol=cv0, dMan=dMan), FUN=sum))
t2 <- reshape(t2, timevar = "dMan", idvar = c("stem", "si", "vol"), direction="wide", sep="")
t1 <- merge(t1, t2, all.x=T)
t1[is.na(t1)] <- 0
t1[,4:NCOL(t1)] <- t1[,4:NCOL(t1)] / t1$all
t1$all <- NULL
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), t1, all.x=T)
#Fill up the missing positions of the matrix
for(i in 4:NCOL(t1)) {
  a <- loess(t1[,i] ~ t1$vol + t1$stem + t1$si, weights=t1$area, control = loess.control(surface = "direct"))
  t1[is.na(t1[,i]),i] <- pmax(0,pmin(1,predict(a, t1)))[is.na(t1[,i])]
}
#Ensure a sum of 1
t1[rowSums(t1[,4:NCOL(t1)]) <= 0, 4:NCOL(t1)] <- 1
t1[,4:NCOL(t1)] <- t1[,4:NCOL(t1)] / rowSums(t1[,4:NCOL(t1)])
#Round
t1[,4:NCOL(t1)] <- t(apply(t1[,4:NCOL(t1)],  1, function(x) roundt(x=x, target=1, digits=2, closest=F, min=0, max=1)))
write.table(t1, "./dat/actprobs.txt")
rm(t1)

#Transition Probabilities: Development without management
#Here you could have the problem that some regions of the grid have low to no
#observation. This could be filled up with estimates from another forest growth
#model
#write.table(with(dat[dat$man=="nomgmt",], data.frame(si=cs, vol0=cv0, vol1=cv1, stem0=cstem0, stem1=cstem1)), "./dat/dNoman.txt", row.names = FALSE)
t1 <- with(dat[dat$man=="nomgmt",], aggregate(list(area=area), list(si=cs, vol0=cv0, vol1=cv1, stem0=cstem0, stem1=cstem1), FUN=sum))
t1$area <- round(t1$area) #round / ceil
t1 <- t1[t1$area > 0,]
t1 <- merge(t1, expand.grid(sapply(nClasses, function(x) 1:x)), all.y=T, by.x=c("si", "vol0", "stem0"), by.y=c("si", "vol", "stem"))
t1$area[is.na(t1$area)] <- 1
t1$vol1[is.na(t1$vol1)] <- t1$vol0[is.na(t1$vol1)]
t1$stem1[is.na(t1$stem1)] <- t1$stem0[is.na(t1$stem1)]
write.table(t1[rep(seq_len(nrow(t1)), t1$area),1:5], "./dat/dNoman.txt", row.names = FALSE)
rm(t1)
#
#estimate transition probabilities
##Create estimation control file
###data ./dat/dNoman.txt .. filename with data from above
###Prior .. could be "uninformative" or a filename wher the matrix is read in
##si=1 .. the order of factors, "from most influential to least influential", 
##together with corresponding weights
cat("data ./dat/dNoman.txt
prior uninformative
si=1
", file="./dat/estiminput.txt")
##levels of the classes
writeLines(sapply(names(nClasses), function(x) {paste(x,paste(1:nClasses[x],collapse = " "))}), "./dat/factors.txt")
##
source("../../EFDMcode/nea/efdmsetuptools.r")
source("../../EFDMcode/nea/efdmestim.r")
#./dat/estiminput.txt .. inputfilename: the name of the estimation control file (that would be given in the activities control file)
#./dat/factors.txt .. factorsfilename: the name of the state space defining control file (SEfactors.txt in the toy dataset)
#vol age .. Changing: names of the factors that transit (given in the Activities control file)
#./dat/pNomgmt.RData .. resultfilename: name of the file where the result is to be saved (do use the extension .RData if you plan to read this in during model/simulation run!)
pre.estimate("./dat/estiminput.txt", "./dat/factors.txt", "vol stem", "./dat/pNomgmt.RData")

#Transition probability for clearcut
#Everything goes back to the first diameter class
t1 <- prod(nClasses[c("vol","stem")])
cat(rep(1, t1), file="./dat/pClearcut.txt")
replicate(t1-1, cat("\n0", rep(0, t1-1), file="./dat/pClearcut.txt", append = T))
cat("\n", file="./dat/pClearcut.txt", append = T)

#Transition probability for man? activities (Default dstem1=destm0, dcv1=dcv0)
for(i in unique(dat$dMan)) {
  shiftarray <- expand.grid(vol=1:nClasses["vol"], stem=1:nClasses["stem"])
  t1 <- unique(dat[dat$dMan==i,c("cv0","cstem0","cv1","cstem1")])
  shiftarray[match(paste(t1$cv0,t1$cstem0), paste(shiftarray$vol,shiftarray$stem)),] <- t1[c("cv1", "cstem1")]
  makeanewP("./dat/pNomgmt.RData", shiftarray, nClasses[c("vol","stem")], paste0("./dat/pMan", i, '.RData'))
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
t1 <- with(dat, aggregate(list(x=drain*area, area=area), list(vol=cv0, stem=cstem0, si=cs), FUN=sum))
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
