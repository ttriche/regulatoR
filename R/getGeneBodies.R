getGeneBodies <- function(SE) { # {{{
  genome <- unique(genome(rowData(SE)))
  if(!(genome %in% c('hg18','hg19'))) {
    stop('Genomes other than hg18 and hg19 are unsupported')
  } else if(genome=='hg19') {
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if(genome=='hg18') {
    require(TxDb.Hsapiens.UCSC.hg18.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
  }
  GR <- reduce(transcriptsBy(txdb, 'gene'))
  SE[ queryHits(findOverlaps(rowData(SE), GR)), ] 
} # }}}
