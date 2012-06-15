alignedToRPKM <- function(readcounts) { # assumes the counts came from Rsubread!

  millionsMapped <- colSums(readcounts$counts)/1000000
  if('ExonLength' %in% names(readcounts$annotation)) {
    geneLengthsInKB <- readcounts$annotation$ExonLength/1000
  } else { 
    geneLengthsInKB <- readcounts$annotation$GeneLength/1000
  }

  # example usage: readcounts$RPKM <- alignedToRPKM(readcounts)
  return( sweep(readcounts$counts, 2, millionsMapped, '/') / geneLengthsInKB )

}
