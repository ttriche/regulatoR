# from exonlvl_exonanaly_hg19.TXT.gz files (ENCODE)
#
# setwd("/where/you/keep/your/processed/exon/files")
exonArrayAsSummarizedExperiment <- function(exlvl.gz) {

  # Windows = crap
  require(parallel)

  # process the ranges specified for the exon array and its probesets
  exlvl = read.delim(exlvl.gz[[1]], stringsAsFactors=F)[, c(1,6,8,9,7)]
  names(exlvl) = c('name','chr','start','end','strand') 
  bogus = which( exlvl$chr %in% c('---','') )
  exlvl = exlvl[ -bogus, ]
  exlvl$strand = as.factor(exlvl$strand)
  levels(exlvl$strand)[ which(levels(exlvl$strand)=='---') ] = '*'
  exlvl.GR = df2GR(exlvl)
  genome(exlvl.GR) = 'hg19'

  # matrix of intensities
  signal <- do.call(cbind, mclapply(exlvl.gz, function(x) {
    read.delim(x, stringsAsFactors=FALSE)[,2]
  }))[ -bogus, ] # drop rows not in the GR
  IDs = sapply(exlvl.gz, function(x) {
    paste(strsplit(x,'_',fixed=T)[[1]][1:2], collapse='_')
  }) # name the columns accordingly
  colnames(signal) = IDs

  EXON.se = SummarizedExperiment(assays=SimpleList(signal=signal), 
                                 colData=DataFrame(sampleNames=IDs),
                                 rowData=exlvl.GR)
  colnames(EXON.se) = EXON.se$sampleNames # don't know why I have to do this...
  rm(exlvl) # the data.frame
  rm(signal) # the matrix
  return(EXON.se)
} # }}}
