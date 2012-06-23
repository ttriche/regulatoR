getRefSeqTxS <- function(SE, what='TES', genome='hg19') { # {{{
  if(!(genome %in% c('hg18','hg19'))) {
    stop('Genomes other than hg18 and hg19 are unsupported')
  } else {
    target <- toupper(paste('refseq', what, genome, sep='.'))
    data(target)
    GR <- get(target)
  }
  SE[ queryHits(findOverlaps(rowData(SE), GR)), ] 
} # }}}

getRefSeqTES <- function(SE, genome='hg19') getRefSeqTxS(SE, 'TES', genome)  
getRefSeqTSS <- function(SE, genome='hg19') getRefSeqTxS(SE, 'TSS', genome)  

getFANTOMTSS <- function(SE, genome='hg19') { # {{{
  if(!(genome %in% c('hg18','hg19'))) {
    stop('Genomes other than hg18 and hg19 are unsupported')
  } else if(genome=='hg18') {
    require(FDb.FANTOM4.promoters.hg19)
    GR <- features(FDb.FANTOM4.promoters.hg19)
  } else if(genome=='hg19') {
    require(FDb.FANTOM4.promoters.hg18)
    GR <- features(FDb.FANTOM4.promoters.hg18)
  }
  SE[ queryHits(findOverlaps(rowData(SE), GR)), ] 
} # }}}
