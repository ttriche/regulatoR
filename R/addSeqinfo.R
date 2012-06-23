addSeqinfo <- function(x) { 
  require(rtracklayer)
  stopifnot( is(x, 'SummarizedExperiment') || is(x, 'GenomicRanges'))
  gen <- ifelse(is(x, 'SummarizedExperiment'), genome(rowData(x)), genome(x))
  gen <- unique(gen)
  if(is.na(gen)) stop('You must assign genome(x) first')
  if(is(x, 'SummarizedExperiment')) {
    seqinfo(rowData(x)) <- SeqinfoForBSGenome(gen)[seqlevels(rowData(x))]
  } else {
    seqinfo(x) <- SeqinfoForBSGenome(gen)[seqlevels(x)]
  }
  return(x)
}
