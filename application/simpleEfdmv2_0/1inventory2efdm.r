#https://ec.europa.eu/jrc/en/european-forestry-dynamics-model/download/required-files

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
funGwl <- function(age, si) { #Estimate total increment until age for si
  c <- 100*si / (150 * si * (1 - exp((-0.028 + 0.00048*si)*100))^(13.5*si^-0.57))
  c * 150 * si * (1 - exp((-0.028 + 0.00048*si)*age))^(13.5*si^-0.57)
}
funStock <- function(age, si, thin=250) {
  funGwl(age, si) / (1 + age/thin)
}
funUopt <- function(age, si) {funGwl(age, si)/age}
funIncRed <- function(age, si, stock) {
  tt <- funStock(age, si)
  bg <- pmin(1, ifelse(tt>0, stock / tt, 1))
  bg^pmax(.2, pmin(1, (1 - 2*si/age)))^2
}
dat <- data.frame(plotId = 1:n, si=round(runif(n,4,16),1), man=factor("nomgmt", levels = c("nomgmt","thin","clearcut")), area=round(runif(n,0.05,10),2))
t1 <- sapply(dat$si, function(x) optimize(funUopt, interval=c(15,199), maximum=T, si=x)$maximum) #Rotation time to reach GWLmax
t1 <- t1 * rnorm(n,mean=1,sd=.2) #Add some nois
dat$age0 <- round(runif(n,0,t1))
dat$age1 <- dat$age0 + dt
dat$vol0 <- funStock(dat$age0, dat$si, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$si/40),.85)
t2 <- (funGwl(dat$age1, dat$si) - funGwl(dat$age0, dat$si)) * funIncRed(dat$age0, dat$si, dat$vol0) #increment
t2 <- t2 * rnorm(n, mean = 1, sd = .1) #Add some nois
dat$vol1 <- dat$vol0 + t2
t3 <- floor(runif(n,0,dt))  #Years until managent action
t4 <- (dat$vol0 + (t2*t3/dt)) - funStock(dat$age0+t3, dat$si, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$si/40),0.85) #Thin more in young wiht high si
dat$man[dat$vol0 + t2 > funStock(dat$age1, dat$si, 1000) * runif(n,0.85,1) & t4 > 0] <- "thin"
dat$man[dat$age1 > t1] <- "clearcut"
dat$drain <- 0
tt <- dat$man=="clearcut"
dat$drain[tt] = dat$drain[tt] + dat$vol0[tt] + (t2*t3/dt)[tt]
dat$age1[tt] <- (dt - t3 - floor(runif(n,0,4)))[tt]
dat$vol1[tt] <- (funStock(pmax(0,dat$age1), dat$si, 1000) * runif(n,0.5,1))[tt]
tt <- dat$man=="thin"
dat$drain[tt] <- dat$drain[tt] + t4[tt]
dat$vol1[tt] <- ((dat$vol0 + (t2*t3/dt)) - t4)[tt]
dat$vol1[tt] <- (dat$vol1 + (funGwl(dat$age1, dat$si) - funGwl(dat$age0+t3, dat$si)) * funIncRed(dat$age0+t3, dat$si, dat$vol1) * rnorm(n, mean = 1, sd = .1))[tt]
dat$vol0 <- round(dat$vol0, 1)
dat$vol1 <- round(dat$vol1, 1)
dat$drain <- round(dat$drain, 1)

#Define ranges of classes
#large volume span to avoid hiting the “roof” of the volume dimension
#probability of growing out of a state for each +- same
#number of volume classes should be equal to the number of time periods (age classes) that it takes to reach its maximum
classVol <- round(funGwl(seq(0, max(c(dat$age0, dat$age1)+dt), dt), weighted.mean(dat$si, dat$area)))
classVol <- round(c(classVol, seq(max(classVol), max(c(dat$vol0, dat$vol1))+2, length.out=2)[-1]))
classAge <- seq(0,max(c(dat$age0,dat$age1))+dt,dt) #Age classes (should coincides with dt), gk:and max should be considered
classSi <- floor(seq(min(dat$si), max(dat$si)+2, 2)) #Siteindex classes
nClasses <- c(length(classVol), length(classAge), length(classSi))-1
names(nClasses) <- c("vol", "age", "si")

#Classify the dataset
dat$ca0 <- findInterval(dat$age0, classAge)
dat$ca1 <- findInterval(pmax(0,dat$age1), classAge)
dat$cs <- findInterval(dat$si, classSi)
dat$cv0 <- findInterval(dat$vol0, classVol)
dat$cv1 <- findInterval(dat$vol1, classVol)

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)


#Initial state from where the simulatin starts: How much area is in the defined classes
#I take here timestep 0, bit also timestep 1 could be used alternatively
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), with(dat, aggregate(list(area=area), list(age=ca0, si=cs, vol=cv0), FUN=sum)), all.x=T)
t1$area[is.na(t1$area)] <- 0
write.table(t1, "./dat/initstate.txt", row.names = FALSE)
rm(t1)

#Activities: transition probabilities depending on different activities which
#are defined here
cat("nomgmt read ../dat/pNomgmt.RData vol age
thin read ../dat/pThin.RData vol age
clearcut read ../dat/pClearcut.txt vol age
", file="./dat/activities.txt")
#no management could be estimated on the fly with "nomgmt estimate estiminput.txt vol age"

#Propabilities for harvest activities
#proportions of the area in a state-space cell that must be subject to the respective activities within a given simulation period
t1 <- with(dat, aggregate(list(all=area, thin=area*(man=="thin"), clearcut=area*(man=="clearcut")), list(age=ca0, si=cs, vol=cv0), FUN=sum))
t1 <- merge(expand.grid(sapply(nClasses, function(x) 1:x)), with(t1, data.frame(age, si, vol, thin=round(thin/all,2), clearcut=round(clearcut/all,2), area=all)), all.x=T)
#Fill up teh matrix with the average, where we have no observations
#Maybe find better solution for that
#t1$clearcut[is.na(t1$clearcut)] <- round(weighted.mean(t1$clearcut, t1$area, na.rm=T),2)
#t1$thin[is.na(t1$thin)] <- round(weighted.mean(t1$thin, t1$area, na.rm=T),2)
#t1$nomgmt <- pmax(0, 1 - (t1$clearcut + t1$thin))
#Fill up teh matrix with loess using the classes
t1$nomgmt <- pmax(0, 1 - (t1$clearcut + t1$thin))
a <- loess(clearcut ~ vol + age + si, data=t1, weights=area, control = loess.control(surface = "direct"))
t1$clearcut[is.na(t1$clearcut)] <- pmax(0,pmin(1,predict(a, t1)))[is.na(t1$clearcut)]
a <- loess(thin ~ vol + age + si, data=t1, weights=area, control = loess.control(surface = "direct"))
t1$thin[is.na(t1$thin)] <- pmax(0,pmin(1,predict(a, t1)))[is.na(t1$thin)]
a <- loess(nomgmt ~ vol + age + si, data=t1, weights=area, control = loess.control(surface = "direct"))
t1$nomgmt[is.na(t1$nomgmt)] <- pmax(0,pmin(1,predict(a, t1)))[is.na(t1$nomgmt)]
#Ensure a sum of 1
t1[rowSums(t1[c("thin", "clearcut", "nomgmt")])<=0,c("thin", "clearcut", "nomgmt")] <- 1
t1[c("thin", "clearcut", "nomgmt")] <- t1[c("thin", "clearcut", "nomgmt")] / rowSums(t1[c("thin", "clearcut", "nomgmt")])
write.table(t1[c("vol", "age", "si", "thin", "clearcut", "nomgmt")], "./dat/actprobs.txt")
rm(t1)

#Transition Probabilities: Development without management
#Here you could have the problem that some regions of the grid have low to no
#observation. This could be filled up with estimates from another forest growth
#model
write.table(with(dat[dat$man=="nomgmt",], data.frame(si=cs, vol0=cv0, vol1=cv1, age0=ca0, age1=ca1)), "./dat/dNoman.txt", row.names = FALSE)
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
source("../../EFDMcode/v2.0/hackfunctions.r")
source("../../EFDMcode/v2.0/efdmestim.r")
#./dat/estiminput.txt .. inputfilename: the name of the estimation control file (that would be given in the activities control file)
#./dat/factors.txt .. factorsfilename: the name of the state space defining control file (SEfactors.txt in the toy dataset)
#vol age .. Changing: names of the factors that transit (given in the Activities control file)
#./dat/pNomgmt.RData .. resultfilename: name of the file where the result is to be saved (do use the extension .RData if you plan to read this in during model/simulation run!)
pre.estimate("./dat/estiminput.txt", "./dat/factors.txt", "vol age", "./dat/pNomgmt.RData")

#Transition probability for thinning activities
#How many volume classes is a forest gowing down when it is thinned
t1 <- with(dat[dat$man=="thin",], aggregate(list(x=(findInterval(vol0, classVol) - findInterval(pmax(0, vol0-drain), classVol))*area, area=area), list(cv=cv0), FUN=sum))
#Thinning volume could be more than vol0
t1$t <- t1$x/t1$area
t1 <- merge(data.frame(cv=1:nClasses["vol"]), t1, all.x=T)
#Fill up mussing data
a <- lm(t ~ cv, data=t1, weights=area)
t1$t[is.na(t1$t)] <- predict(a, t1)[is.na(t1$t)]
t1$t[t1$t>t1$cv-1] <- t1$cv[t1$t>t1$cv-1]-1 #Ensure drops are in range
t1$t[t1$t<0] <- 0 #No negativ numbers
#./dat/pNomgmt.RData .. inputfilename: name of the save file where the R object for "no management" matrix is
#t1 .. voldrops: a vector of integers, how many levels to drop at each level (length equal to number of volume levels)
#./dat/pThin.RData .. resultfilename: name of the file where the result is to be saved (do use the extension .RDat
#makeathinP("./dat/pNomgmt.RData", round(t1$t), "./dat/pThin.RData") #It looks like that it works also with float numbers
makeathinP("./dat/pNomgmt.RData", t1$t, "./dat/pThin.RData")

#Transition probability for clearcut
#Everything goes back to the first age class
t1 <- prod(nClasses[c("vol","age")])
cat(rep(1, t1), file="./dat/pClearcut.txt")
replicate(t1-1, cat("\n0", rep(0, t1-1), file="./dat/pClearcut.txt", append = T))
cat("\n", file="./dat/pClearcut.txt", append = T)



#Output request file
cat("volume by age
../dat/orVolumeAge.txt
drain by age si
../dat/orDrainAgeSi.txt
area by age
../dat/orAreaAge.txt
area by age si
../dat/orAreaAgeSi.txt
drain by vol age si
../dat/orDrainAgeSi.txt
volume by vol age si
../dat/orVolumeVolAgeSi.txt
", file="./dat/outputrequests.txt")
#Average Volume per volume and age class
t1 <- with(dat, aggregate(list(x=c(vol0,vol1)*area, area=area), list(vol=c(cv0,cv1), age=c(ca0,ca1)), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + age, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","age")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, age, volume=round(t,1))), "./dat/orVolumeAge.txt", row.names = FALSE)
#Average Harvest per volume, age and si class
t1 <- with(dat, aggregate(list(x=drain*area, area=area), list(vol=cv0, age=ca0, si=cs), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + age + si, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","age","si")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, age, si, drain=round(t,1))), "./dat/orDrainAgeSi.txt", row.names = FALSE)
#Area per age
t1 <- data.frame(age = 1:nClasses["age"], area=1)
write.table(t1, "./dat/orAreaAge.txt", row.names = FALSE)
#Area per age and si class
t1 <- expand.grid(sapply(nClasses[c("age","si")], function(x) 1:x))
t1$area <- 1
write.table(t1, "./dat/orAreaAgeSi.txt", row.names = FALSE)
#Average Volume per volume, age and si class
t1 <- with(dat, aggregate(list(x=c(vol0,vol1)*area, area=area), list(vol=c(cv0,cv1), age=c(ca0,ca1), si=c(cs,cs)), FUN=sum))
t1$t <- t1$x/t1$area
a <- loess(t ~ vol + age + si, data=t1, control = loess.control(surface = "direct"))
t1 <- merge(expand.grid(sapply(nClasses[c("vol","age","si")], function(x) 1:x)), t1, all.x=T)
t1$t[is.na(t1$t)] <- pmin(max(t1$t, na.rm=T), pmax(0, predict(a, newdata=t1[is.na(t1$t),])))
write.table(with(t1, data.frame(vol, age, si, volume=round(t,1))), "./dat/orVolumeVolAgeSi.txt", row.names = FALSE)


#Create efdm control file
cat("STATESPACE_FILENAME ../dat/factors.txt
INITSTATE_FILENAME ../dat/initstate.txt
ACTIVITIES_FILENAME ../dat/activities.txt
ACTPROBS_FILENAME ../dat/actprobs.txt
NROFSTEPS 3
OUTPUTREQUEST_FILENAME ../dat/outputrequests.txt
", file="./dat/efdminput.txt")
