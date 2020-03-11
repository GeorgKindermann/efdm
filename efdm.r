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
  if(sparse) {   #Sparse array
    x <- aggregate(x[,i]
           , list(i=efdmIndexGet(efdmClassGet(x[!i], TRUE, breaks), t0)), sum)
    x <- Matrix::sparseVector(x$x, x$i, prod(t0))
  } else {     #Dense array
    x <- as.vector(xtabs(reformulate(termlabels = ".", response = what), data.frame(x[i], efdmClassGet(x[!i], FALSE, breaks))))
  }
  attr(x, "dimS") <- t0
  x
}

efdmTransitionGet <- function(x, area="area", flow=c("harvest","residuals"), t0="0", t1="1", breaks=breaks, model=c("qda","lda"), dimWgt=1) {
  env <- new.env(parent=emptyenv())
  env$dimS <- efdmDimGet(colnames(x)[endsWith(colnames(x), t0)], breaks)
  env$dimWgt <- rep_len(dimWgt, length(env$dimS))
  x <- aggregate(. ~ from + to, data.frame(from = efdmIndexGet(efdmClassGet(x[,endsWith(colnames(x), t0)], TRUE, breaks), env$dimS)
   , to = efdmIndexGet(efdmClassGet(x[,endsWith(colnames(x), t1)], TRUE, breaks), env$dimS)
   , tmp[match(c(area,flow), colnames(tmp))]), sum)
  tmp <- aggregate(area~from, x, sum)
  x$share <- x$area / tmp$area[match(x$from, tmp$from)]
  x[flow] <- x[flow] / x[,area]
  env$simple <- x[c("from","to","share",flow)]
  if(!is.null(model) && !is.na(model[1])) {
    fun <- function(i, j) {
      match.fun(model[j])(reformulate(termlabels = names(env$dimS)
                , response = paste0(i,".1")), tt, weights= x[area])}
    tt <- data.frame(efdmIndexInv(x$from,env$dimS),efdmIndexInv(x$to,env$dimS))
    env$model <- lapply(setNames(names(env$dimS), names(env$dimS)), function(i){
      for(j in seq_along(model)) {
        try(ret <- fun(i,j), TRUE)
        if(exists("ret", inherits = FALSE)) return(ret)
      }})
  }
  env
}

efdmNextArea <- function(t0, transition, share="share") {
  t1 <- if(is.vector(t0)) {numeric(length(t0))
        } else {Matrix::sparseVector(0, 1, length(t0))}
  attr(t1, "dimS") <- attr(t0, "dimS")
  static <- efdmIndexStatic(setdiff(names(attr(t0, "dimS")), names(transition$dimS)), attr(t0, "dimS"))
  from <- efdmIndexGet(efdmIndexInv(transition$simple$from, transition$dimS), attr(t0, "dimS"))
  to <- efdmIndexGet(efdmIndexInv(transition$simple$to, transition$dimS), attr(t0, "dimS"))
  tt <- aggregate(as.vector(t0[rowSums(expand.grid(from, static))]) * transition$simple[,share], list(to=rowSums(expand.grid(to, static))), sum)
  tt <- tt[tt$x > 0,]
  t1[tt$to] <- t1[tt$to] + tt$x
  i <- if(is.vector(t0)) {which(t0 != 0)} else {t0@i}
  i <- i[!i %in% rowSums(expand.grid(from, static))]
  if(length(i) > 0) {
    if(exists("model", envir = transition)) {
      from <- as.data.frame(efdmIndexInv(i, attr(t0, "dimS")))
      to <- do.call(cbind, lapply(transition$model, function(f) {
        x <- predict(f, from)$class #Maybe move to more than one cell
        as.numeric(levels(x))[x] }))
      to <- efdmIndexGet(cbind(to, from[setdiff(colnames(from), colnames(to))])
                       , attr(t0, "dimS"))
      from <- efdmIndexGet(from, attr(t0, "dimS"))
      tt <- aggregate(as.vector(t0[from]), list(to=to), sum)
      t1[tt$to] <- t1[tt$to] + tt$x
    } else {t1[i] <- t1[i] + t0[i]}
  }
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
