keepFeature <- function(x, ..., lg=FALSE) { 
  # this takes less than 0.08 seconds/feature on my laptop 
  rows <- rownames(x) 
  if(is(x, 'SummarizedExperiment')) {
    x <- asy.fast(x)
    rownames(x) <- rows
  }
  if(lg == TRUE) x <- log1p(x)
  stopifnot(is(x, 'matrix'))
  require(parallel)
  unlist(mclapply(rows, function(y) return(mixFilter(x[y,],...)$retain > 1)))
}
