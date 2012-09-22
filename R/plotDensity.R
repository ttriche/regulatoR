plotDensity <- function(x, minSep=0.2, ..., lg=F) { 
  rows <- rownames(x)
  if(is(x, 'SummarizedExperiment')) x <- asy.fast(x)
  if(lg) x <- log1p(x)
  rownames(x) <- rows
  x <- x[ !is.na(x) ]
  res = clusterAndCollapse(x, minSep=minSep, ...)
  means = unique(res$values)
  plot(density(x, adjust=0.5), , lty=1, lwd=2, ..., 
       main=paste(res$retain,'components (dashed) at minSep =',minSep))
  if(res$retain > 0) {
    for(i in seq_along(means)) abline(v=means[i], col=i+2, lty=3, lwd=2)
    dens=density(res$values, adjust=0.5)
    lines(dens$x, dens$y*(1.2/length(means)), lty=3, lwd=2)
  }
}
