transMatWP <- function(x, weight, prior) {
  freq <- xtabs(area~., x, addNA=T)
  t1 <- sapply(x[,-NCOL(x)], nlevels)
  tt <- c(t1[1]*t1[2], t1[3]*t1[4])
  if(length(t1>4)) {tt <- c(tt, t1[5:length(t1)])}
  dim(freq) <- tt
  s <- prior
  for(i in 2:length(tt)) {
    tmp <- apply(freq, 1:i, sum)
    #if(weight[i-1] == 1) {
      t1 <- which(colSums(tmp) == 0, arr.ind = TRUE)
      t1 <- matrix(t1, ncol=NCOL(t1))
      tmp[cbind(rep(seq_len(NROW(tmp)), each = NROW(t1)), t1[rep(1:nrow(t1), NROW(tmp)), ])] <- s[cbind(rep(seq_len(NROW(tmp)), each = NROW(t1)), t1[rep(1:nrow(t1), NROW(tmp)), seq_len(max(1,NCOL(t1)-1))])]
    #}
    tmp <- prop.table(tmp,2:length(dim(tmp)))
    #tmp[is.na(tmp)] <- 0
    s <- sweep(tmp*weight[i-1], 1:max(2,i-1), s*(1-weight[i-1]), '+')
    #s <- prop.table(s,2:length(dim(s)))
  }
  #prop.table(s,2:length(tt))
  s
}

getPriorIJ <- function(i, j, di, dj) {
  ret <- matrix(0, nrow=i*j, ncol=i*j)
  tt <- expand.grid(i=1:i,j=1:j)
  tt$ti <- pmax(1, pmin(i, tt$i+di))
  tt$tj <- pmax(1, pmin(j, tt$j+dj))
  ret[cbind((tt$tj-1)*i+tt$ti, (tt$j-1)*i+tt$i)] <- 1
  ret
}

getPriorObs <- function(x, di=0, dj=0) {
  tt <- expand.grid(levels(x[,3]), levels(x[,4]))
  names(tt) <- names(x[,3:4])
  me <- merge(x, tt, all.y=T)
  me <- me[names(x)]
  me[is.na(me[,NCOL(me)]), NCOL(me)] <- 1
  
  tt[,names(x[,1:2])] <- tt
  tt[,3] <- levels(tt[,3])[pmax(1, pmin(nlevels(tt[,3]), match(tt[,3], levels(tt[,3])) + di))]
  tt[,4] <- levels(tt[,4])[pmax(1, pmin(nlevels(tt[,4]), match(tt[,4], levels(tt[,4])) + di))]
  tt <- tt[names(x)[1:4]]
  t1 <- !complete.cases(me[,1:2])
  me[t1,1:2] <- tt[match(interaction(me[t1,3:4]), interaction(tt[,3:4])),1:2]

  freq <- xtabs(as.formula(paste0(names(me)[NCOL(me)], "~.")), me[,c(1:4,NCOL(me))], addNA=T)
  dim(freq) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]))
  prop.table(freq,2)
}

require(randomForest)
getPriorRF <- function(x) {
  tt <- data.frame(x[,1:2], apply(x[,c(3:4,NCOL(x))], 2, as.numeric))
  a1 <- randomForest(as.formula(paste0(names(tt)[1], " ~ ", names(tt)[3], " + ", names(tt)[4])), data=tt, weights=tt$area)
  a2 <- randomForest(as.formula(paste0(names(tt)[2], " ~ ", names(tt)[3], " + ", names(tt)[4])), data=tt, weights=tt$area)

  tt <- expand.grid(as.numeric(levels(x[,3])), as.numeric(levels(x[,4])))
  names(tt) <- names(x[,3:4])

  t1 <- predict(a1, tt, type = "prob")
  t2 <- predict(a2, tt, type = "prob")
  sapply(seq_len(NROW(t1)), function(i) {t1[i,] * rep(t2[i,], each=NCOL(t1))})
}

require(e1071)
getPriorSvm <- function(x) {
  tt <- data.frame(x[,1:2], apply(x[,c(3:4,NCOL(x))], 2, as.numeric))
  a1 <- svm(as.formula(paste0(names(tt)[1], " ~ ", names(tt)[3], " + ", names(tt)[4])), data=tt, weights=tt$area, probability = T)
  a2 <- svm(as.formula(paste0(names(tt)[2], " ~ ", names(tt)[3], " + ", names(tt)[4])), data=tt, weights=tt$area, probability = T)

  tt <- expand.grid(as.numeric(levels(x[,3])), as.numeric(levels(x[,4])))
  names(tt) <- names(x[,3:4])

  t1 <- attr(predict(a1, tt, probability = T), "probabilities")
  t2 <- attr(predict(a2, tt, probability = T), "probabilities")
  sapply(seq_len(NROW(t1)), function(i) {t1[i,] * rep(t2[i,], each=NCOL(t1))})
}

#transMatRF <- function(x) {
#}

#transMatSvm <- function(x) {
#}
