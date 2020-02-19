transMatWP <- function(x, weight, prior) {
  freq <- xtabs(area~., x, addNA=T)
  t1 <- sapply(x[,-NCOL(x)], nlevels)
  tt <- c(t1[1]*t1[2], t1[3]*t1[4])
  if(length(t1)>4) {tt <- c(tt, t1[5:length(t1)])}
  dim(freq) <- tt
  s <- prior
  for(i in length(dim(prior)):length(tt)) {
    tmp <- apply(freq, 1:i, sum)
    #if(weight[i-1] == 1) {
      t1 <- which(colSums(tmp) == 0, arr.ind = TRUE)
      if(NROW(t1) > 0) {
      t1 <- matrix(t1, ncol=NCOL(t1))
      tmp[cbind(rep(seq_len(NROW(tmp)), each = NROW(t1)), t1[rep(1:nrow(t1), NROW(tmp)), ])] <- s[cbind(rep(seq_len(NROW(s)), each = NROW(t1)), t1[rep(1:nrow(t1), NROW(s)), seq_len(min(length(dim(s))-1,NCOL(t1)))])]
      }
    #}
    tmp <- prop.table(tmp,2:length(dim(tmp)))
    #tmp[is.na(tmp)] <- 0
    s <- sweep(tmp*weight[i-1], 1:length(dim(s)), s*(1-weight[i-1]), '+')
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
  tt <- expand.grid(sapply(3:(NCOL(x)-1), function(i) levels(x[,i])))
  names(tt) <- names(x[,3:(NCOL(x)-1)])
  me <- merge(x, tt, all.y=T)
  me <- me[names(x)]
  me[is.na(me[,NCOL(me)]), NCOL(me)] <- 1

  tt[,names(x[,1:2])] <- tt[,1:2]
  tt <- tt[names(x)[1:(NCOL(x)-1)]]
  tt[,1] <- levels(tt[,1])[pmax(1, pmin(nlevels(tt[,1]), match(tt[,1], levels(tt[,1])) + di))]
  tt[,2] <- levels(tt[,2])[pmax(1, pmin(nlevels(tt[,2]), match(tt[,2], levels(tt[,2])) + dj))]
  t1 <- !complete.cases(me[,1:2])
  me[t1,1:2] <- tt[match(interaction(me[t1,3:4]), interaction(tt[,3:4])),1:2]
  
  freq <- xtabs(as.formula(paste0(names(me)[NCOL(me)], "~.")), me, addNA=T)
  if(NCOL(x)>5) {
    dim(freq) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]), sapply(5:(NCOL(x)-1), function(i) nlevels(x[,i])))
  } else {
    dim(freq) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]))
  }
  prop.table(freq,2:length(dim(freq)))
}

if(require(MASS)) {
getPriorLda <- function(x) {
  tt <- data.frame(x[1:2], do.call(cbind, lapply(x[,c(3:NCOL(x))], function(x) unclass(x))))
  a1 <- lda(as.formula(paste0(names(tt)[1], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)])
  a2 <- lda(as.formula(paste0(names(tt)[2], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)])
  tt <- expand.grid(sapply(3:(NCOL(x)-1), function(i)  seq_along(levels(x[,i]))))
  names(tt) <- names(x[,3:(NCOL(x)-1)])
  t1 <- predict(a1, tt)$posterior
  t2 <- predict(a2, tt)$posterior
  tt <- sapply(seq_len(NROW(t1)), function(i) {t1[i,] * rep(t2[i,], each=NCOL(t1))})
  if(NCOL(x)>5) {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]), sapply(5:(NCOL(x)-1), function(i) nlevels(x[,i])))
  } else {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]))
  }
  tt
}
}

if(require(randomForest)) {
getPriorRF <- function(x) {
  tt <- data.frame(x[1:2], apply(x[,c(3:NCOL(x))], 2, as.numeric))
  a1 <- randomForest(as.formula(paste0(names(tt)[1], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)])
  a2 <- randomForest(as.formula(paste0(names(tt)[2], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)])

  tt <- expand.grid(sapply(3:(NCOL(x)-1), function(i) as.numeric(levels(x[,i]))))
  names(tt) <- names(x[,3:(NCOL(x)-1)])

  t1 <- predict(a1, tt, type = "prob")
  t2 <- predict(a2, tt, type = "prob")
  tt <- sapply(seq_len(NROW(t1)), function(i) {t1[i,] * rep(t2[i,], each=NCOL(t1))})
  if(NCOL(x)>5) {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]), sapply(5:(NCOL(x)-1), function(i) nlevels(x[,i])))
  } else {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]))
  }
  tt
}
}

if(require(e1071)) {
getPriorSvm <- function(x) {
  tt <- data.frame(x[1:2], apply(x[,c(3:NCOL(x))], 2, as.numeric))
  a1 <- svm(as.formula(paste0(names(tt)[1], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)], probability = T)
  a2 <- svm(as.formula(paste0(names(tt)[2], " ~ ", paste(names(tt)[3:(NCOL(tt)-1)], collapse=" + "))), data=tt, weights=tt[NCOL(tt)], probability = T)

  tt <- expand.grid(sapply(3:(NCOL(x)-1), function(i) as.numeric(levels(x[,i]))))
  names(tt) <- names(x[,3:(NCOL(x)-1)])

  t1 <- attr(predict(a1, tt, probability = T), "probabilities")
  t2 <- attr(predict(a2, tt, probability = T), "probabilities")
  tt <- sapply(seq_len(NROW(t1)), function(i) {t1[i,] * rep(t2[i,], each=NCOL(t1))})
  if(NCOL(x)>5) {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]), sapply(5:(NCOL(x)-1), function(i) nlevels(x[,i])))
  } else {
    dim(tt) <- c(nlevels(x[,1])*nlevels(x[,2]), nlevels(x[,3])*nlevels(x[,4]))
  }
  tt
}
}

transMat2LT <- function(x, mat) {
  tt <- data.frame(
    factor(levels(x[,1])[1 + (slice.index(mat, 1)-1) %% nlevels(x[,1])], levels=levels(x[,1]), ordered=is.ordered(x[,1])),
    factor(levels(x[,2])[1 + floor((slice.index(mat, 1)-1) / nlevels(x[,1]))], levels=levels(x[,2]), ordered=is.ordered(x[,2])),
    factor(levels(x[,3])[1 + (slice.index(mat, 2)-1) %% nlevels(x[,3])], levels=levels(x[,3]), ordered=is.ordered(x[,3])),
    factor(levels(x[,4])[1 + floor((slice.index(mat, 2)-1) / nlevels(x[,3]))], levels=levels(x[,4]), ordered=is.ordered(x[,4]))
  )
  if(NCOL(x)>5) {
    tt <- cbind(tt, as.data.frame(lapply(5:(NCOL(x)-1), function(i) factor(levels(x[,i])[slice.index(mat, i-2)], levels=levels(x[,i]), ordered=is.ordered(x[,i])))))
  }
  tt <- cbind(tt, as.vector(mat))
  names(tt) <- names(x)
  tt <- tt[tt[,NCOL(tt)] > 0,]
  tt
}

if(require(MASS)) {
  fillEmptyTarget <- function(x) {
        res <- x[0,]
    for(j in 1:2) {
      t1 <- as.data.frame(setdiff(levels(x[,j]), x[,j]))
      if(NROW(t1) > 0) {
        colnames(t1) <- names(x)[j]
        t2 <- x[0,][seq_len(NROW(t1)),]
        t2[,j] <- factor(t1[,1], levels=levels(t2[,j]), ordered=is.ordered(t2[,j]))
        t2[,NCOL(t2)] <- mean(x[,NCOL(x)])/NROW(t2)/100
        for(i in setdiff(seq_len(NCOL(x)-1), j)) {
          t2[,i] <- tryCatch({
            a <- lda(as.formula(paste0(names(x)[i], " ~ as.numeric(", names(x)[j],")")), data=x)
            predict(a, t1)$class
          }, error = function(error_condition) {
            names(which.max(table(x[,names(x)[i]])))
          })
        }
        res <- rbind(res, t2)
      }
    }
    res
  }
}

getTrans <- function(x, tre=0) {
  sapply(1:2, function(i) {
    t1 <- sort(xtabs(as.formula(paste0(names(x)[NCOL(x)], " ~ I(match(", names(x)[i], ", levels(", names(x)[i], ")) - match(", names(x)[i+2], ", levels(", names(x)[i+2], ")))")), data=x))
    as.integer(names(t1)[cumsum(t1) > sum(t1*tre)])
  }, simplify = F)
}

cutLT <- function(lt, trans) {
  res <- lt
  for(i in 1:2) {
    tt <- unique(data.frame(levels(lt[,i+2]), levels(lt[,i+2])[pmax(1, pmin(nlevels(lt[,i+2]), seq_len(nlevels(lt[,i+2])) + rep(trans[[i]], each=nlevels(lt[,i+2]))))]))
    names(tt) <- names(lt)[c(i+2,i)]
    res <- merge(res, tt)
  }
  res[names(lt)]
}

#library(mda)
#mda(vol~.,data=x)
#fda(vol~.,data=x)

#library(MASS)
#qda(vol~.,data=x)

#library(klaR)
#rda(vol~.,data=x)

#library(nnet)
#nnet(vol~.,data=x)

#library(kernlab)
#ksvm(vol~.,data=x)

#library(caret)
#knn3(vol~.,data=x)

#library(e1071)
#naiveBayes(vol~.,data=x)

