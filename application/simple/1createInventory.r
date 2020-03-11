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
