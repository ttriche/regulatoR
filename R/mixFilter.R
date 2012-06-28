mixFilter <- function(x, minSep=0.2, G=c(2:5), ...) {
  bak <- list(values=x, retain=0)
  res <- try(clusterAndCollapse(x, minSep=minSep, G=G), silent=TRUE)
  if(!inherits(res, 'try-error')) {
    return(res)
  } else {
    return(bak)
  }
}
