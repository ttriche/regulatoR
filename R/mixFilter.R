mixFilter <- function(x, ...) {
  bak <- list(values=x, retain=0)
  res <- try(clusterAndCollapse(x, ...), silent=T)
  if(!inherits(res, 'try-error')) return(res)
  else return(bak)
}
