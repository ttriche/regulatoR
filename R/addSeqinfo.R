addSeqinfo <- function(SE) { 
  require(rtracklayer)
  if(any(is.na(unique(genome(SE))))) stop('You must assign a genome to your SE')
  seqinfo(SE) <- SeqinfoForBSGenome(unique(genome(SE)))[seqlevels(SE)]
  SE
}
