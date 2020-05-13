#Schaftholz nach Braun: (d>=10.6cm)
fzBraunCl <- function(coefBraun) {
  force(coefBraun)
  fun <- function(typ, a, b1, b2, b3, d, h, d03h, hk) {
    if(is.na(typ) | !(typ > 0 & typ < 5)) {return(NA)}
    switch(typ
         , conif={a + b1*d03h/d + b2*h/d + b3/d}
         , fir={a + b1*(d03h/d)^2 + b2*h/d + b3/d}
         , broadl={a + b1*d03h/d + b2*hk/h + b3*h/d^2}
         , oak={a + b1*(d03h/d)^2 + b2*hk/h + b3*h/d^2})
  }
  function(x) {
    unlist(with(cbind(x[,c("d","h","d03h","hk")], coefBraun[match(x$bart, coefBraun$s),c("t","a","b1","b2","b3")]), mapply(fun,t,a,b1,b2,b3,d/10,h*10,d03h/10,hk*10)))
  }
}
#
#Koeffizienten Braun für Messungen in dm
#wird in Funktion auf d[cm], d03h[cm], h[m], hk[m] umgerechnet
coefBraun <- read.table(header = TRUE, text = "
t s   a        b1       b2          b3
1  1 -0.243560 0.827140 0.000291030 0.0287120   #Fichte
2  2  0.099066 0.512570 0.000445934 0.0159819   #Tanne
1  3 -0.219830 0.802760 0.000323901 0.0183610   #Laerche
1  4 -0.209860 0.814030 0.000195795 0.0317280   #Kiefer
1  5 -0.192940 0.847910 0.000203583 0.0069137   #Schw.Kiefer
1  6  0.050144 0.467565 0.000157070 0.0761060   #Zirbe
3 10 -0.130890 0.674300 0.066763    0.000167307 #Buche
4 11  0.185205 0.350116 0.065650    0.000476876 #Eiche
3 12  0.042132 0.422578 0.077020    0.000420750 #Hainbuche
3 13 -0.019811 0.512416 0.053539    0.000470395 #Esche
3 14 -0.028594 0.565491 0.008268    0.000236824 #Ahorn
3 15 -0.138981 0.695015 0.016567    0.000318013 #Ulme
3 20 -0.077804 0.568230 0.051675    0.000554142 #Birke
3 21 -0.164590 0.703820 0.058915    0.000259377 #Erle
3 23 -0.143768 0.648730 0.081205    0.000561830 #Weisspappel
3 24 -0.084312 0.592805 0.022711    0.000647000 #Schwarzpappel
3 25 -0.145570 0.665700 0.058867    0.000417790 #Zitter+Baösampappel
3 27 -0.137580 0.694380 0.012841    0.000458940 #Weide
")

#Schaftholz nach Schieler (5 <= d <10.6cm)
fzSchielCl <- function(coefSchiel) {
  force(coefSchiel)
  fun <- function(b1, b2, b3, b4, b5, b6, b7, d, h) {
    b1 + b2*log(d)^2 + b3/h + b4/d + b5/d^2 + b6/(d*h) + b7/(d^2*h)
  }
  function(x) {
    unlist(with(cbind(x[,c("d","h")], coefSchiel[match(x$bart, coefSchiel$s), c("b1","b2","b3","b4","b5","b6","b7")]), mapply(fun,b1,b2,b3,b4,b5,b6,b7,d/10,h*10)))
  }
}
#
#Koeffizienten Schieler für Messungen in dm 
#wird in Funktion auf d[cm], d03h[cm], h[m], hk[m] umgerechnet
coefSchiel <- read.table(header = TRUE, text = "
s b1        b2        b3       b4        b5       b6       b7
 1 0.563443 -0.12731  -8.55022  0        0        7.6331 0       #Fichte
 2 0.560673  0.15468  -0.65583  0.03321  0        0      0       #Tanne
 3 0.48727   0        -2.04291  0        0        5.9995 0       #Laerche
 4 0.435949 -0.0149083 5.21091  0        0.028702 0      0       #Weisskiefer
 5 0.53438  -0.00763   0        0        0        0      2.2414  #Schwarzkiefer
 6 0.525744 -0.0334896 7.38943 -0.10646  0        0      3.34479 #Zirbe
10 0.5173    0       -13.62144  0        0        9.9888 0       #Buche
11 0.417118  0.21941  13.32594  0        0        0      0       #Eiche
12 0.32473   0.02432   0        0.23972  0       -9.9388 0       #Hainbuche
13 0.48122  -0.01489 -10.83056  0        0        9.3936 0       #Esche
14 0.50101  -0.03521  -8.07176  0        0.03521  0      0       #Ahorn
15 0.44215  -0.02446   0        0        0        0      2.87714 #Ulme
20 0.42831  -0.06643   0        0        0        8.4307 0       #Birke
21 0.387399  0         7.17123  0.04407  0        0      0       #Erle
24 0.366419  0         1.13323  0.1306   0        0      0       #Pappel
25 0.4115   -0.00989 -28.27478  0.35599 -0.21986 21.4913 0       #Schwarzpappel
27 0.54008  -0.02716 -25.11447  0.08327  0        9.3988 0       #Weide
")

bam <- read.csv("./dat/bam.csv.bz2")
names(bam)[match(c("bhd", "hoehe", "kronho"), names(bam))] <- c("d", "h", "hk")
bam$d <- bam$d/10
bam$d03h <- bam$d03h/10
bam$h <- bam$h/10
bam$hk <- bam$hk/10
bam$nrep <- 16/(bam$d/100)^2/pi
bam$nrep[bam$d < 10.4] <- 10000/(2.6^2*pi)

tt <- data.frame(s=sort(unique(bam$bart)))
tt$braun <- coefBraun$s[match(tt$s, coefBraun$s)]
tt$schiel <- coefSchiel$s[match(tt$s, coefSchiel$s)]

tt$schiel[tt$s == 23] <- 21
tt[tt$s %in% c(7.0, 8.0, 9.0, 9.1, 9.2, 9.3, 9.5, 9.6), c("braun", "schiel")] <- c(4,1,1,1,2,3,4,2)
tt[tt$s %in% c(16.0, 17.0, 18.0, 18.1, 18.2, 18.4, 22.0, 24.1, 24.2, 26.0, 31.0, 31.1, 31.2, 31.4, 32.0), c("braun", "schiel")] <- c(11,11,10,10,10,10,21,24,24,24,11,11,11,11,27)

library(compiler)
#fzBraun <- cmpfun(fzBraunCl(merge(tt[,c("s","braun")], coefPoll, by.x="braun", by.y="s")[-1]))
fzBraun <- cmpfun(fzBraunCl(cbind(data.frame(s=tt$s), subset(coefBraun, select=-c(s))[match(tt$braun, coefBraun$s),])))
#fzSchiel <- cmpfun(fzBraunCl(merge(tt[,c("s","schiel")], coefPoll, by.x="schiel", by.y="s")[-1]))
fzSchiel <- cmpfun(fzSchielCl(cbind(data.frame(s=tt$s), subset(coefSchiel, select=-c(s))[match(tt$schiel, coefSchiel$s),])))

bam$fz <- fzBraun(bam)
key <- is.finite(bam$d) & is.finite(bam$h) & (bam$d<10.4 | is.na(bam$fz))
bam$fz[key] <- fzSchiel(bam[key, c("bart", "d", "h")])

#stem(bam$fz)
#quantile(bam$fz[bam$peri==6], na.rm=T, probs = seq(0, 1, 0.025))
#with(bam[bam$peri==6,], summary(do.call("rbind", tapply(fz, bart, quantile, na.rm=T, probs = c(0.025, 0.975)))))
bam$fz[bam$fz < .25] <- .25
bam$fz[bam$fz > .75] <- .75
bam$fz[is.na(bam$fz)] <- with(bam[bam$peri==6,], mean(fz, na.rm=T, weight=nrep*d^2*h*fz))
bam$v <- with(bam, (d/200)^2*pi*h*fz*nrep)
bam$v[is.na(bam$v)] <- 0
key <- c("nutzart","ea", "ba")
bam[key] <- lapply(bam[key], factor)

prf <- read.csv("./dat/prf.csv.bz2")
prf$datum <- as.Date(prf$datum)

tfl <- read.csv("./dat/tfl.csv.bz2")
key <- c("ba","ea","relief","wugeb","bogrup","vegtyp")
tfl[key] <- lapply(tfl[key], factor)

si <- read.csv("./dat/oberhoeheFichte.csv.bz2")

t1 <- prf[prf$peri %in% 6:7, c("rw","hw","pbfl","peri","datum")]
tt <- as.numeric(t1$datum - as.Date(paste0(format(t1$datum,"%Y"),"-05-01"))) / 122
tt[tt < 0] <- 0
tt[tt > 1] <- 1
t1$vergYear <- as.numeric(format(t1$datum,"%Y")) + tt
t1 <- reshape(t1, timevar="peri", idvar = c("rw","hw","pbfl"), direction="wide", drop="datum")
summary(t1$vergYear.7 - t1$vergYear.6)
#Im Mittel 7 Jahre

tt <- bam[,c("peri","rw","hw","pbfl","tlfl","azi","dist","bart","nutzart","v","d")]
tt <- reshape(tt, timevar="peri", idvar = c("rw","hw","pbfl","tlfl","azi","dist"), direction="wide")
gc()
tt$bart <- ifelse(is.na(tt$bart.7), tt$bart.6, tt$bart.7)
tt$nutzart <- tt$nutzart.7
tt$d0 <- tt$d.6
tt$d1 <- tt$d.7
tt$v0 <- tt$v.6
tt$v1 <- tt$v.7
tt <- tt[!(is.na(tt$d0) & is.na(tt$d1)),!grepl("\\.[67]", colnames(tt))]

tt <- merge(prf[prf$peri==6 & prf$kg_wald>0, c("rw", "hw", "pbfl", "kg_wald")], tt, all.x=TRUE)
tt$v0[is.na(tt$v0)] <- 0
tt$v1[is.na(tt$v1)] <- 0
tt[rowSums(is.na(tt[,c("d0", "d1")])) == 2, c("d0", "d1")] <- 0

tt$nutzart[is.na(tt$nutzart)] <- ""
table(tt$nutzart[!is.na(tt$d1)])
tt$nutzart[tt$nutzart=="" & !is.na(tt$d0) & is.na(tt$d1)] <- NA

t1 <- aggregate(cbind(v = v0+v1) ~ rw + hw + bart, tt[!is.na(tt$bart),], sum)
t1 <- t1[order(t1$rw, t1$hw, -t1$v),]
t1 <- aggregate(bart ~ rw + hw, t1, "[", 1)
i <- is.na(tt$bart)
j <- match(interaction(tt[i,c("rw","hw")]), interaction(t1[,c("rw","hw")]))
tt$bart[i] <- t1$bart[match(interaction(tt[i,c("rw","hw")]), interaction(t1[,c("rw","hw")]))]
i <- which(is.na(tt$bart))
tt$bart[i] <- sapply(i, function(i) t1$bart[which.min(abs(tt$rw[i] - t1$rw) + abs(tt$hw[i] - t1$hw))])

tt <- merge(tt, si, all.x=TRUE)
t3 <- aggregate(si ~ rw + hw, data=si, FUN=mean)
t2 <- is.na(tt$si)
tt$si[t2] <- t3$si[match(interaction(tt[t2,c("rw","hw")]), interaction(t3[,c("rw","hw")]))]
a <- loess(si ~ rw + hw, data=si, control=loess.control(surface="direct"))
t2 <- is.na(tt$si)
tt$si[t2] <- predict(a, newdata=tt[t2,c("rw","hw")])
rm(a)

t1 <- with(prf, data.frame(rw,hw,pbfl,seehoehe,rueck,vorrueck,peri,dis=abs(peri-6)))
t1 <- t1[order(t1$dis, -t1$peri),]
t1 <- aggregate(cbind(seehoehe, rueck, vorrueck) ~ rw + hw + pbfl, t1, function(x) x[which.min(is.na(x))], na.action=NULL)
tt <- merge(tt, t1, all.x=TRUE)

t1 <- unique(tt[,c("rw", "hw", "pbfl", "tlfl")])
t1 <- merge(t1, tfl[tfl$peri==6, c("rw","hw","pbfl","tlfl","ba","ea","hangneig","neigricht","relief","wugeb","bogrup","vegtyp")], all.x=TRUE)
for(i in c("ba","ea","hangneig","neigricht","relief","wugeb","bogrup","vegtyp")) {
  for(peri in c(6:max(tfl$peri), 5:min(tfl$peri))) {
    j <- is.na(t1[,i])
    k <- match(interaction(t1[j,c("rw","hw","pbfl")]), interaction(tfl[tfl$peri==peri, c("rw","hw","pbfl")]))
    t1[j,i] <- tfl[tfl$peri==6, i][k]
  }
}
tt <- merge(tt, t1, all.x=TRUE)

t1 <- lapply(tt[, (match("si", names(tt))+1):ncol(tt)], function(x) {
  if(any(is.na(x))) {
    i <- is.na(x)
    if(is.factor(x)) {x <- as.character(x)}
    x[i] <- type.convert(names(sort(table(x[!i]), decreasing = TRUE))[1], as.is = TRUE)
    x
  } else {x}
  })
tt <- cbind(tt[, 1:match("si", names(tt))], as.data.frame(t1))

t1 <- aggregate(cbind(ant0=v0, ant1=v1) ~ rw + hw + pbfl, data=tt, FUN=sum)
tt <- merge(tt, t1, all.x=TRUE)
tt$ant0 <- tt$v0 / tt$ant0
tt$ant1 <- tt$v1 / tt$ant1
tt$ant0[is.na(tt$ant0)] <- 1
tt$ant1[is.na(tt$ant1)] <- 1
tt$ant0[is.na(tt$d0)]  <- 0
tt$ant1[is.na(tt$d1)]  <- 0
t1 <- do.call(data.frame, aggregate(cbind(ant0, ant1) ~ rw + hw + pbfl, tt, function(x) c(sum=sum(x), n=length(x))))
tt <- merge(tt, t1[round(t1$ant0.sum) == 0, c("rw","hw","pbfl","ant0.n")], all.x=TRUE)
tt <- merge(tt, t1[round(t1$ant1.sum) == 0, c("rw","hw","pbfl","ant1.n")], all.x=TRUE)
i <- !is.na(tt$ant0.n) & is.na(tt$ant1.n)
tt$ant0[i] <- tt$ant1[i]
tt$d0[i] <- 0
i <- !is.na(tt$ant1.n) & is.na(tt$ant0.n)
tt$ant1[i] <- tt$ant0[i]
tt$d1[i] <- 0
#i <- !is.na(tt$ant0.n) & !is.na(tt$ant1.n) #Maybe nott needed
#tt$ant0[i] <- 1/tt$ant0.n[i]
#tt$ant1[i] <- 1/tt$ant1.n[i]
tt$ant0.n <- NULL
tt$ant1.n <- NULL

tt$area <- 55*27.5/4 * tt$kg_wald / 10

sum(tt$area * tt$ant0, na.rm=TRUE)
sum(tt$area * tt$ant1, na.rm=TRUE)

write.table(tt[,c("rw","hw","pbfl","bart","nutzart","d0","d1","v0","v1","ant0","ant1","area","si","seehoehe","rueck","vorrueck","ba","ea","hangneig","neigricht","relief","wugeb","bogrup","vegtyp")], file=bzfile("/tmp/owiEfdmBase.txt.bz2"), row.names=FALSE)


#mv /tmp/owiEfdmBase.txt.bz2 ./dat/

