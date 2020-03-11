x <- read.csv("./dat/inventoryData.csv.bz2")
y <- read.table("./res/rawoutput.txt")

t1 <- with(x, aggregate(list(obs0=area), list(vol=vol0, stem=stem0, si=si, man=man), FUN=sum))
t1 <- reshape(t1, v.names = "obs0", timevar = "man", idvar = c("vol", "stem", "si"), direction = "wide")
t2 <- with(x, aggregate(list(obs1=area), list(vol=vol1, stem=stem1, si=si), FUN=sum))

me <- merge(merge(t1, t2, all=T), y, all=T)
me[is.na(me)] <- 0
me$est1 <- rowSums(me[,grep("step1", colnames(me))])

colSums(me[,4:NCOL(me)])

sunflowerplot(me$obs0.nomgmt, me$step0.nomgmt)
sunflowerplot(me$obs1, me$est1)

with(me, boxplot(obs1 - est1 ~ vol))
with(me, boxplot(obs1 - est1 ~ stem))
with(me, boxplot(obs1 - est1 ~ si))

t1 <- with(me, aggregate(list(obs1=obs1, nomgmt.1=est1), list(vol=vol, stem=stem), FUN=sum))

tt <- matrix(0, nrow = max(me$stem), ncol = max(me$vol))
tt[cbind(t1$stem, t1$vol)] <- t1$obs1
image(1:max(me$stem), 1:max(me$vol), tt, xlab="Stem Class", ylab="Volume Class", main="Observed Area")
with(t1, text(stem, vol, round(obs1)))

tt <- matrix(0, nrow = max(me$stem), ncol = max(me$vol))
tt[cbind(t1$stem, t1$vol)] <- t1$nomgmt.1
image(1:max(me$stem), 1:max(me$vol), tt, xlab="Stem Class", ylab="Volume Class", main="Predicted Area")
with(t1, text(stem, vol, round(nomgmt.1,2)))

tt <- matrix(0, nrow = max(me$stem), ncol = max(me$vol))
tt[cbind(t1$stem, t1$vol)] <- t1$obs1 - t1$nomgmt.1
image(1:max(me$stem), 1:max(me$vol), tt, xlab="Stem Class", ylab="Volume Class", main="Observed - Estimated Area")
with(t1, text(stem, vol, round(obs1 - nomgmt.1,2)))

