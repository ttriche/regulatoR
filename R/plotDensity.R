plotDensity <-
function(x, minSep=0.1, xlab='Value', more='', ...) { 
  x <- x[ !is.na(x) ]
  res = clusterAndCollapse(x, minSep)
  means = unique(res$values)
  plot(density(x, adjust=0.5), xlab=xlab, lty=1, lwd=2, ..., 
       main=paste0(res$retain,' components (dashed) at minSep=',minSep,more))
  if(res$retain > 0) {
    for(i in seq_along(means)) abline(v=means[i], col=i+2, lty=3, lwd=2)
    dens=density(res$values, adjust=0.5)
    lines(dens$x, dens$y*(1.2/length(means)), lty=3, lwd=2)
  }
}
