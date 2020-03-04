require(Matrix)
#library(slam) #Alternatives for spares arrays
#library(spray)

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
