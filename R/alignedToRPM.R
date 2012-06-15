alignedToRPM <- function(readcounts) { # assumes the counts came from Rsubread!
  millionsMapped <- colSums(readcounts$counts)/1000000
  return( sweep(readcounts$counts, 2, millionsMapped, '/') )
}
