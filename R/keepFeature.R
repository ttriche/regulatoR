keepFeature <- function(x, minSep=0.2, G=c(1:5), ...) {
  stopifnot(is(x, 'matrix'))
  require(parallel)
  unlist(mclapply(rownames(x), function(y) {
    return( mixFilter(x[y,],minSep=minSep,G=G)$retain > 1 )
  }))
}
