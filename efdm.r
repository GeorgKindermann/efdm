LIBdistances <- require(distances)
LIBslam <- require(slam)

require(Matrix)
require(MASS)

c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

efdmBreakNamesGet <- function(x, env=breaks) {
  i <- do.call(rbind, lapply(ls(breaks), function(y) {
    i <- grep(paste0("^", y), x)
    if(length(i)>0) {data.frame(i, y, stringsAsFactors = FALSE)}}))
  if(length(x) == nrow(i)) {i[order(i[,1]),]}
}

efdmClassGet <- function(x, int=TRUE, env=breaks) {
  i <- efdmBreakNamesGet(colnames(x), env)
  as.data.frame(lapply(seq_len(nrow(i)), function(j) {
    k <- findInterval(x[,j], get(i[j,2], env))
    if(int) {k
    } else {factor(k,0:length(get(i[j,2], env)))}
  }), col.names=colnames(x))
}

efdmDimGet <- function(x, env=breaks) {
  vapply(unique(efdmBreakNamesGet(x, env)[,2]), function(i) {length(get(i, env))}, numeric(1))+1
}

efdmIndexGet <- function(x, d) {
  i <- lapply(paste0("^", names(d)), grep, colnames(x))
  as.vector(matrix(unlist(x), nrow(x)) %*% cumprod(c(1,d))[match(seq_len(ncol(x)), i)]) + 1
}

efdmIndexInv <- function(i, d) {
  m <- unname(rev(cumprod(d))[-1])
  "colnames<-"(do.call(rbind, lapply(i, function(i) rev(Reduce("%%", m, (i-1), accumulate = TRUE) %/% c(m, 1)))), names(d))
}

efdmIndexStatic <- function(i, d) {
  if(length(i) == 0) {0
  } else {
    if(is.character(i)) {i <- match(i, names(d))}
    m <- c(1, cumprod(d))
    rowSums(expand.grid(lapply(i, function(i) seq(0, by=m[i], length=d[i]))))
  }
}

efdmStateGet <- function(x, what="area", breaks=breaks, sparse=TRUE) {
  i <- colnames(x) == what
  t0 <- efdmDimGet(colnames(x)[!i], breaks)
  if(sparse && LIBslam) {   #Sparse array
    x <- aggregate(x[,i]
           , list(i=efdmIndexGet(efdmClassGet(x[!i], TRUE, breaks), t0)), sum)
    x <- slam::simple_sparse_array(x$i, x$x, prod(t0))
  } else {     #Dense array
    x <- as.vector(xtabs(reformulate(termlabels = ".", response = what), data.frame(x[i], efdmClassGet(x[!i], FALSE, breaks))))
  }
  attr(x, "dimS") <- t0
  x
}

efdmTransitionGet <- function(x, area="area", t0="0", t1="1", breaks=breaks, getArea=TRUE) {
  env <- new.env(parent=emptyenv())
  c0 <- grep(paste0("^", area, "|", t1, "$"), colnames(x), value = TRUE, invert =TRUE)
  c1 <- grep(paste0("^", area, "|", t0, "$"), colnames(x), value = TRUE, invert =TRUE)
  env$dimS <- efdmDimGet(c0, breaks)
  x <- aggregate(list(area=x[,area]), list(from = efdmIndexGet(efdmClassGet(x[c0], TRUE, breaks), env$dimS), to = efdmIndexGet(efdmClassGet(x[c1], TRUE, breaks), env$dimS)), sum)
  tmp <- aggregate(area~from, x, sum)
  x$share <- x$area / tmp$area[match(x$from, tmp$from)]
  env$simple <- cbind(x[c("from","to","share")], area= x[,area[getArea]])
  env
}

efdmTransitionXda <- function(env, model=c("qda","lda"), absolute=FALSE, dropArea=TRUE) {
  fun <- function(i, j) {
    match.fun(model[j])(reformulate(termlabels = names(env$dimS)
    , response = paste0(i,".1")), tt, weights= env$simple$area)}
  tt <- data.frame(efdmIndexInv(env$simple$from,env$dimS), efdmIndexInv(env$simple$to,env$dimS))
  if(!absolute) {
    i <- colnames(tt) %in% names(env$dimS)
    tt[!i] <- tt[!i] - tt[i]
  }
  env$model <- lapply(setNames(names(env$dimS), names(env$dimS)), function(i){
    for(j in seq_along(model)) {
      try(ret <- fun(i,j), TRUE)
      if(exists("ret", inherits = FALSE)) return(ret)
    }})
  env$absolute <- absolute
  if(dropArea) {env$simple$area <- NULL}
  env
}

efdmTransitionSimilar <- function(env, absolute=TRUE, wgt=1, fixed=FALSE) {
  env$absolute <- absolute
  env$dimWgt <- rep_len(wgt, length(env$dimS))
  env$similar <- rep_len(!fixed, length(env$dimS))
}

efdmTransitionShifter <- function(env, shift=0) {
  env$shifter <- rep_len(shift, length(env$dimS))
}


efdmNextArea <- function(t0, transition, nMinToGo=1, nMaxToGo=3, pMinToGo=0.05, maxNeighbours=7) {
  staticN <- setdiff(names(attr(t0, "dimS")), names(transition$dimS))
  static <- efdmIndexStatic(staticN, attr(t0, "dimS"))
  from <- efdmIndexGet(efdmIndexInv(transition$simple$from, transition$dimS), attr(t0, "dimS"))
  to <- efdmIndexGet(efdmIndexInv(transition$simple$to, transition$dimS), attr(t0, "dimS"))
  res <- aggregate(list(area=t0[rowSums(expand.grid(from, static))] * transition$simple[,"share"]), list(from=rowSums(expand.grid(from, static)), to=rowSums(expand.grid(to, static))), sum) #Suboptimal
  rm(from, to)
  res <- res[res$area > 0,]
  i <- if(is.vector(t0)) {which(t0 != 0)} else {t0$i[,1]}
  i <- i[!i %in% res$from]
  if(length(i) > 0) {
    if(exists("model", envir = transition)) {
      from <- as.data.frame(efdmIndexInv(i, attr(t0, "dimS")))
      to <- lapply(setNames(seq_along(transition$model), names(transition$model)), function(i) {
        f <- transition$model[[i]]
        if(is.null(f)) {setNames(rep(1,nrow(from)),from[,names(transition$model)[i]])
        } else {
          x <- predict(f, from)$posterior
          rownames(x) <- NULL
          if(ncol(x) > 1) {
            apply(x, 1, function(x) proportions(head(x[order(-x)], min(nMaxToGo, max(nMinToGo, sum(x >= pMinToGo))))))
          } else {asplit(x, 1)}
        }})
      x <- do.call(cbind, to)
      to <- do.call(cbind, lapply(to, function(x) if(is.list(x)) lapply(x, function(x) as.integer(names(x))) else as.integer(names(x))))
      if(!transition$absolute) {
        to <- do.call(cbind, lapply(setNames(colnames(to), colnames(to)), function(j) {tt <- unlist(Map("+", from[,j], to[,j]))
          k <- tt >= transition$dimS[j]
          tt[k] <- (transition$dimS - 1)[j]
          k <- tt < 0
          tt[k] <- 0
          relist(tt, to[,j])}))
      }
      rm(from)
      x <- apply(x, 1, function(x) apply(expand.grid(x), 1, prod))
      to <- apply(to, 1, function(x) efdmIndexGet(expand.grid(x), attr(t0, "dimS")))
      x <- mapply("*", x, t0[i])
      res <- rbind(res, aggregate(list(area=unlist(x)), list(from=rep(i, lengths(to)), to=unlist(to)), sum))
    } else if(exists("similar", envir = transition)) {
      from <- as.matrix(as.data.frame(efdmIndexInv(i, attr(t0, "dimS"))))
      from <- from[,names(transition$dimS)]
      fromObsI <- unique(transition$simple$from)
      fromObs <- as.matrix(as.data.frame(efdmIndexInv(fromObsI, transition$dimS)))
      x <- rbind(fromObs, from)
      x <- distances(x * transition$dimWgt[col(x)])
      idx <- seq_len(nrow(from))+nrow(fromObs)
      to <- nearest_neighbor_search(x, maxNeighbours, idx, seq_len(length(fromObsI)))
      to <- lapply(seq_along(idx), function(j) {
        y <- x[idx[j],to[,j]]
        to[,j][y == y[1]]
      })
      rm(x, idx)
      transition$simple <- transition$simple[order(transition$simple$from),]
      x <- cbind(1L+findInterval(fromObsI[unlist(to)]-1, transition$simple$from), findInterval(fromObsI[unlist(to)], transition$simple$from))
      x <- mapply(seq, x[,1], x[,2])
      to <- unname(split(x, rep(seq_along(to), lengths(to))))
      rm(x)
      to <- lapply(to, function(y) {j <- unlist(y)
        list(transition$simple$to[j], proportions(transition$simple$share[j])
           , transition$simple$from[j]) #Only needed for relativ shift
      })
      to <- do.call(rbind, to)
      colnames(to) <- c("to", "share", "from")
      to[,"share"] <- mapply("*", to[,"share"], t0[i])
      if(exists("absolute", envir = transition) & transition$absolute) {
        to[,"to"] <- mapply(function(x,y) {
          tt <- efdmIndexInv(x, transition$dimS)
          tt[!transition$similar[col(tt)]] <- rep(y[!transition$similar], each=nrow(tt))
          tt[,staticN] <- y[staticN]
          efdmIndexGet(tt, attr(t0, "dimS"))
        }, to[,"to"], split(from, seq_len(nrow(from))))
      } else {
        to[,"to"] <- mapply(function(x,y) {
          tt <- x + y[col(x)]
          j <- tt >= transition$dimS[col(tt)]
          tt[j] <- (transition$dimS - 1)[col(tt)][j]
          j <- tt < 0
          tt[j] <- 0
          tt[,staticN] <- y[staticN]
          efdmIndexGet(tt, attr(t0, "dimS"))}, 
          mapply(function(x,y) {
            efdmIndexInv(x, transition$dimS) - efdmIndexInv(y, transition$dimS)},
            to[,"to"], to[,"from"]), split(from, seq_len(nrow(from))))
      }
      res <- rbind(res, aggregate(list(area=unlist(to[,"share"])), list(from=rep(i, lengths(to[,"to"])), to=unlist(to[,"to"])), sum))
    } else if(exists("shifter", envir = transition)) {
      from <- as.matrix(as.data.frame(efdmIndexInv(i, attr(t0, "dimS"))))
      from <- from[,names(transition$dimS)]
      to <- from + transition$shifter[col(from)]
      j <- to >= transition$dimS[col(to)]
      to[j] <- (transition$dimS - 1)[col(to)][j]
      j <- to < 0
      to[j] <- 0
      to <- efdmIndexGet(to, attr(t0, "dimS"))
      res <- rbind(res, aggregate(list(area=t0[i]), list(from=i, to=to), sum))
    } else {
      res <- rbind(res, data.frame(from=i, to=i, area=t0[i]))
    }
  }
  aggregate(area ~ from + to, res, sum)
}

efdmNextState <- function(t0, x) {
  t1 <- if(is.vector(t0)) {numeric(length(t0))
        } else {simple_sparse_array(1, 0, dim(t0))}
  tt <- aggregate(area ~ to, x, sum)
  t1[tt$to] <- tt$area
  attr(t1, "dimS") <- attr(t0, "dimS")
  t1
}


efdmNextFlow <- function(t0, transition, flow="harvest", share="share") {
  res <- if(is.vector(t0)) {numeric(length(t0))
         } else {Matrix::sparseVector(0, 1, length(t0))}
  attr(res, "dimS") <- attr(t0, "dimS")
  static <- efdmIndexStatic(setdiff(names(attr(t0, "dimS")), names(transition$dimS)), attr(t0, "dimS"))
  from <- efdmIndexGet(efdmIndexInv(transition$simple$from, transition$dimS), attr(t0, "dimS"))
  tt <- aggregate(as.vector(t0[rowSums(expand.grid(from, static))]) * transition$simple[,share] * transition$simple[,flow], list(to=rowSums(expand.grid(from, static))), sum)
  tt <- tt[tt$x > 0,]
  res[tt$to] <- tt$x
  i <- if(is.vector(t0)) {which(t0 != 0)} else {t0@i}
  i <- i[!i %in% rowSums(expand.grid(from, static))]
  if(length(i) > 0) {
    tt <- aggregate(cbind(x = transition$simple[,flow] * transition$simple[,share]), list(from=transition$simple$from), sum)
    tFrom <- efdmIndexInv(tt$from, transition$dimS)
    iFrom <- efdmIndexInv(i, attr(t0, "dimS"))
    iFrom <- iFrom[,match(colnames(tFrom), colnames(iFrom))]
    res[i] <- res[i] + t0[i] * vapply(seq_along(i), function(k) {
      j <- abs(sweep(tFrom, 2, iFrom[k,])) %*% transition$dimWgt
      mean(tt$x[which(j == min(j))])}, numeric(1))
  }
  res
}
