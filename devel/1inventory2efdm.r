set.seed(0)
dir.create("dat")

#Create a simple example dataset
n <- 1000 #Number of forest stands
dt <- 10 #Timespann between two observations
#plotId..Plot identification
#age0..Age at timestep 0
#age1..Age at timestep 1 (typical age0+dt when there was not clearcut
#volume0..Volume [m3/ha] at timestep 0
#volume1..Volume [m3/ha] at timestep 1
#diameter0..Diamter [cm] at timestep 0
#diameter1..Diamter [cm] at timestep 1
#shareC0..area share of coniferous trees [1] at timestep 0
#shareC1..area share of coniferous trees [1] at timestep 1
#yieldC..Yield level coniferous as average total increment until age 100 [m3/ha/year]
#yieldNC..Yield level nonconiferous as average total increment until age 100 [m3/ha/year]
#man..Managementtype: nomgmt,thin,clearcut,selectiveHarvest
#harvest..Harvests between age0 and age1 [m3/ha/dt]
#residuals..Harvest residuals and mortality between age0 and age1 [m3/ha/dt]
#temperature..Average temperature [degree Celsius]
#area..Representig area of the plot [hectar]

#Helper functions
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
    qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
} #Truncated normal distribution https://stackoverflow.com/a/19350640/10488504
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

dat <- data.frame(
  plotId = 1:n
, man=factor("nomgmt", levels =c("nomgmt","thin","clearcut","selectiveHarvest"))
, area = round(runif(n,0.05,10),2)
, shareC0 = round(runif(n,0,1),2)
, temperature = round(runif(n,2,16),1)
)
dat$shareC1 <- pmax(0, round(pmin(1, round(dat$shareC + runif(n,-.2,.2), 2)),2))
dat$yieldC <- pmax(1, round(14 - abs(8 - dat$temperature) + runif(n,-2,2), 1))
dat$yieldNC <- pmax(1, round(14 - abs(12 - dat$temperature) + runif(n,-2,2), 1))
dat$yield <- round(dat$yieldC*dat$shareC0 + dat$yieldNC*(1-dat$shareC0), 1)
t1 <- sapply(dat$yield, function(x) optimize(funUopt, interval=c(15,199), maximum=T, si=x)$maximum) #Rotation time to reach GWLmax
t1 <- t1 * rtnorm(n,mean=1,sd=.2,.3,1.7) #Add some nois
dat$age0 <- round(runif(n,0,t1))
dat$age1 <- dat$age0 + dt
dat$volume0 <- funStock(dat$age0, dat$yield, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$yield/40),.85)
dat$diameter0 <- funD(dat$age0, dat$yield) * rtnorm(n, 1, .1, .5, 1.5)
dat$diameter1 <- dat$diameter0 + (funD(dat$age1, dat$yield) - funD(dat$age0, dat$yield)) * funDAdjBg(dat$age0, dat$yield, funGetBg(dat$age0, dat$yield, dat$volume0), dat$diameter0) * rtnorm(n, 1, .1, .5, 1.5)
t2 <- (funGwl(dat$age1, dat$yield) - funGwl(dat$age0, dat$yield)) * funIncRed(dat$age0, dat$yield, funGetBg(dat$age0, dat$yield, dat$volume0)) #increment
t2 <- t2 * rtnorm(n, mean = 1, sd = .1, .5, 1.5) #Add some nois
dat$volume1 <- dat$volume0 + t2
t3 <- floor(runif(n,0,dt))  #Years until managent action
t4 <- (dat$volume0 + (t2*t3/dt)) - funStock(dat$age0+t3, dat$yield, 1000) * runif(n,pmin(.75, 0.5 + dat$age0/dat$yield/40),0.85) #Thin more in young with high si
dat$man[dat$volume0 + t2 > funStock(dat$age1, dat$yield, 1000) * runif(n,0.85,1) & t4 > 0] <- "thin"
dat$man[dat$age1 > t1] <- sample(c("clearcut", "nomgmt"), size=sum(dat$age1 > t1), replace=T, prob=c(.7,.3))
dat$man[t4 > 0 & runif(n) < .3 & funGetBg(dat$age0, dat$yield, dat$volume0) > .8] <- "selectiveHarvest"
dat$harvest <- 0
dat$residuals <- 0
t5 <- dat$volume0 + (t2*t3/dt)
i <- dat$man=="clearcut"
dat$harvest[i] <- dat$harvest[i] + t5[i]*.7
dat$residuals[i] <- dat$residuals[i] + t5[i]*.3
dat$age1[i] <- pmax(1, round(dt - t3 - floor(runif(n,0,4)))[i])
dat$volume1[i] <- (funStock(pmax(0,dat$age1), dat$yield, 1000) * runif(n,0.5,1))[i]
dat$diameter1[i] <- (funD(dat$age1, dat$yield) * rtnorm(n, 1, .1, .5, 1.5))[i]
i <- dat$man%in%c("thin","selectiveHarvest")
dat$harvest[i] <- dat$harvest[i] + t4[i]*.6
dat$residuals[i] <- dat$residuals[i] + t4[i]*.4
dat$volume1[i] <- ((dat$volume0 + (t2*t3/dt)) - t4)[i]
dat$volume1[i] <- (dat$volume1 + (funGwl(dat$age1, dat$yield) - funGwl(dat$age0+t3, dat$yield)) * funIncRed(dat$age0+t3, dat$yield, funGetBg(dat$age0+t3, dat$yield, dat$volume1)) * rtnorm(n, mean = 1, sd = .1, .5, 1.5))[i]
i <- dat$man=="thin"
dat$diameter1[i] <- pmax(1, dat$diameter1 + rtnorm(n, mean = 2, sd = 1, -2, 4))[i]
i <- dat$man=="selectiveHarvest"
dat$diameter1[i] <- pmax(1, dat$diameter1 + rtnorm(n, mean = -4, sd = 2, -6, 2))[i]
dat$age1[i] <- pmax(1, round(dat$age1 + rtnorm(n, mean = -15, sd = 10, -30, 0))[i])
i <- grep("^(temperature|diameter)", colnames(dat))
dat[,i] <- round(dat[,i],1)
i <- grep("^(volume|harvest|residuals)", colnames(dat))
dat[,i] <- round(dat[,i])
dat$yield <- NULL

#Save the Inventory data
write.csv(dat, bzfile("./dat/inventoryData.csv.bz2"), row.names = FALSE)

#Define ranges of classes
#The size of a growth classes (volumen, age, diameter) should be maximum the change in one growht period, so that all leave the class it in one period
#age
breaks <- new.env(parent=emptyenv())
breaks$age <- seq(10,max(c(dat$age0,dat$age1)),dt)
breaks$volume <- seq(30,max(c(dat$volume0,dat$volume1)),30)
breaks$diameter <- seq(3,max(c(dat$diameter0,dat$diameter1)),3)
breaks$yield <- seq(ceiling(min(dat$yieldC,dat$yieldNC)),max(dat$yieldC,dat$yieldNC),1)
breaks$temperature <- seq(ceiling(min(dat$temperature)),max(dat$temperature),1)
save(breaks, file="./dat/breaks.RData")

#Initial state from where the simulation starts: How much area is in the defined classes
#I take here timestep 0, but also timestep 1 could be used alternatively
#Maximum size .Machine$integer.max
t0 <- tt <- c("volume","diameter","temperature")
t0[1:2]  <- paste0(t0[1:2], 0)
key <- as.data.frame(lapply(setNames(seq_along(tt),tt), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- rep("yield", 2)
t0 <- paste0(tt, c("C", "NC"))
keyY <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
stateInit <- xtabs(area ~ ., cbind(key, yield=c(keyY$yieldC, keyY$yieldNC), species=rep(0:1, each=n), area=c(dat[,"area"] * dat$shareC0, dat[,"area"] * (1-dat$shareC0))))
attr(stateInit, "call") <- NULL
#object.size(stateInit)
save(stateInit, file="./dat/initstate.RData", compress="bzip2")

#Transition Probabilities
t0 <- tt <- c(rep(c("volume","diameter"),2))
t0[1:2]  <- paste0(t0[1:2], 0)
t0[3:4]  <- paste0(t0[3:4], 1)
key <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tt <- rep("yield", 2)
t0 <- paste0(tt, c("C", "NC"))
keyY <- as.data.frame(lapply(setNames(seq_along(t0),t0), function(i) {
  factor(findInterval(dat[,t0[i]], get(tt[i], breaks)), 0:length(get(tt[i], breaks)))
}))
tmp <- aggregate(cbind(areaC0=shareC0*area, areaC1=shareC1*area, area) ~ ., cbind(key, keyY, dat[c("shareC0","shareC1","area")]), sum)
tmp <- as.data.frame(lapply(tmp, function(x) {if(is.factor(x)) as.integer(levels(x))[x] else x}))
t0 <- pmin(1, tmp$areaC1/tmp$areaC0)
t0[is.na(t0)] <- 0
t1 <- pmin(1, (tmp$area-tmp$areaC1)/(tmp$area-tmp$areaC0))
t1[is.na(t1)] <- 0
tmp  <- rbind(data.frame(tmp[1:4], species0 = 0, species1 = 0, yield0=tmp$yieldC, yield1=tmp$yieldC, area = tmp$areaC0 * t0)
, data.frame(tmp[1:4], species0 = 0, species1 = 1, yield0=tmp$yieldC, yield1=tmp$yieldNC, area = tmp$areaC0 * (1-t0))
, data.frame(tmp[1:4], species0 = 1, species1 = 1, yield0=tmp$yieldNC, yield1=tmp$yieldNC, area = (tmp$area - tmp$areaC0) * t1)
, data.frame(tmp[1:4], species0 = 1, species1 = 0, yield0=tmp$yieldNC, yield1=tmp$yieldC, area = (tmp$area - tmp$areaC0) * (1-t1)))
tmp <- aggregate(area~.,tmp[tmp$area>0,],sum)
tt <- lapply(paste0("^", names(dimnames(stateInit))), grep, colnames(tmp))
for(i in seq_len(length(tt))) {
  tmp[,tt[[i]]] <- tmp[,tt[[i]]] * cumprod(c(1,dim(stateInit)))[i]
}
tmp <- data.frame(from = rowSums(tmp[,grep("0$", colnames(tmp))])+1
, to = rowSums(tmp[,grep("1$", colnames(tmp))])+1
, area=tmp$area)
tt <- aggregate(area~from, tmp, sum)
tmp$share <- tmp$area / tt$area[findInterval(tmp$from, tt$from)]
tmp$area <- NULL
tmp$from <- as.integer(tmp$from)
tmp$to <- as.integer(tmp$to)
dynamic <- do.call(rbind, unname(lapply(split(tmp, tmp$to), function(x) data.frame(to=x$to[1], fs=I(list(list(x$from, x$share)))))))
#Alternative structure from tmp maybe smaller but probably slower
#object.size(tmp)
#object.size(dynamic)
#dynamic <- tmp
#t1 <- stateInit #Hot ho come with this from t0 to t1
#t1[rowSums(expand.grid(dynamic$from, static))] <- 0
#for(i in seq_len(nrow(dynamic))) {
#  t1[dynamic$to[i] + static] <- t1[dynamic$to[i] + static] + stateInit[dynamic$from[i] + static] * dynamic$share[i]
#}
rm(tmp)
i <- match(setdiff(names(dimnames(stateInit)), c("volume", "diameter", "species", "yield")), names(dimnames(stateInit)))
static <- as.integer(rowSums(expand.grid(lapply(i, function(i) seq(0, by=cumprod(c(1,dim(stateInit)))[i], length.out=dim(stateInit)[i])))))
transition <- new.env(parent=emptyenv())
transition$dynamic <- dynamic
transition$static <- static
save(transition, file="./dat/transition.RData", compress="xz")

#Go to the next timestep
t0 <- stateInit
t1 <- t0
t1[rowSums(expand.grid(unlist(lapply(dynamic$fs, "[[", 1)), static))] <- 0
t1[rowSums(expand.grid(dynamic$to, static))] <- t1[rowSums(expand.grid(dynamic$to, static))] + unlist(lapply(dynamic$fs, function(i) colSums(matrix(t0[rep(static, each=length(i[[1]])) + i[[1]]] * i[[2]], length(i[[1]])))))

