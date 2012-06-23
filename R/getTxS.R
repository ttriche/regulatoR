getRefSeqTxS <- function(SE, what='TES') { # {{{
  genome <- unique(genome(rowData(SE)))
  if(!(genome %in% c('hg18','hg19'))) {
    stop('Genomes other than hg18 and hg19 are unsupported')
  } else {
    target <- toupper(paste('refseq', what, genome, sep='.'))
    data(target)
    GR <- get(target)
  }
  if(what == 'TES') GR <- resize(GR, 200, fix='center')
  if(what == 'TSS') GR <- resize(resize(GR, 200, fix='end'), 250, fix='start')
  SE[ queryHits(findOverlaps(rowData(SE), GR)), ] 
} # }}}

getRefSeqTES <- function(SE) getRefSeqTxS(SE, what='TES')
getRefSeqTSS <- function(SE) getRefSeqTxS(SE, what='TSS')

getFANTOMTSS <- function(SE) { # {{{
  genome <- unique(genome(rowData(SE)))
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
