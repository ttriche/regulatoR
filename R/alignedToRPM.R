alignedToRPM <- function(readcounts) { # assumes the counts came from Rsubread!
  millionsMapped <- colSums(readcounts$counts)/1000000
  return( sweep(readcounts$counts, 2, millionsMapped, '/') )
}

asRPM <- function(SE) { # for SEs where we have counts and/or RPKM but not RPM
  stopifnot('counts' %in% names(assays(SE, withDimnames=F)))
  millionsMapped <- colSums(assays(SE, withDimnames=F)[['counts']])/1000000
  return( sweep(assays(SE, withDimnames=F)[['counts']], 2, millionsMapped, '/'))
}
