# processing (here I am assuming TCGA patient IDs as names)
#
# for i in *exon*; do 
#   j=`echo $i | cut -f1,3 -d'-' | tr - _`
#   cat $i | cut -f1,4 | gzip > $j.rpkm.gz
# done
#
## FIXME: tabulate raw counts as well: FIXED! 2/22
#
# setwd("/where/you/keep/your/processed/RPKM/files")
exonsAsSummarizedExperiment <- function(exons.gz=NULL, genome.aligned=NULL) {

  if(is.null(exons.gz)) { # {{{ print usage tips 
    message("\n")
    message('How to use exonsAsSummarizedExperiment:')
    message("")
    message('1) process TCGA exon RPKMs with, say, bash:')
    message('$  for i in *exon*; do')
    message(">    j=`echo $i | cut -f1,3 -d'-' | tr - _`")
    message(">    cat $i | cut -f1,3,4 | gzip > $j.rpkm.gz")
    message('>  done')
    message("")
    message("2) read them into a list of files in R:")
    message("R> exons.gz = list.files(pattern='.rpkm.gz')")
    message("")
    message("3) run the function using this list:")
    message("R> exons.RPKM = exonsAsSummarizedExperiment(exons.gz)")
    return(FALSE)
  } # }}}

  # Windows = crap
  require(parallel)

  # process the ranges specified in the file (and GAF)
  exons = read.delim(exons.gz[[1]], stringsAsFactors=F)[,1]
  exons = t(sapply(exons, function(x) strsplit(x, ':', fixed=TRUE)[[1]]))
  exons = cbind(exons[,1], 
                t(sapply(exons[,2], function(x) strsplit(x,'-',fixed=T)[[1]])),
                exons[,3])
  rownames(exons) = 1:nrow(exons)
  exons = as.data.frame(exons)
  exons[,2:3] = apply(exons[,2:3], 2, as.numeric)
  names(exons) = c('chr','start','end','strand')

  # matrix of counts
  counts <- do.call(cbind, mclapply(list.files(patt='rpkm.gz$'), function(x) {
    read.delim(x, stringsAsFactors=FALSE)[,2]
  }))
  # matrix of RPKM
  RPKM <- do.call(cbind, mclapply(list.files(patt='rpkm.gz$'), function(x) {
    read.delim(x, stringsAsFactors=FALSE)[,3]
  }))
  # vector of sampleNames
  IDs = unlist(lapply(list.files(patt='rpkm.gz$'), function(x) {
    strsplit(x, '.', fixed=T)[[1]][1]
  }))
  colnames(RPKM) = colnames(counts) = IDs
  EXONS.se = SummarizedExperiment(assays=SimpleList(RPKM=RPKM, counts=counts), 
                                  colData=DataFrame(sampleNames=IDs),
                                  rowData=df2GR(exons))
  colnames(EXONS.se) = EXONS.se$sampleNames # I don't know why I have to do this
  rm(exons) # the data.frame
  rm(counts) # the matrix
  rm(RPKM) # the matrix
  if(is.null(genome.aligned)) {
    message('Be sure to set genome(rowData(your.exons)) and assign $gene_id!')
  } else { 
    genome(rowData(EXONS.se)) <- genome.aligned
  }
  return(EXONS.se[ order(rowData(EXONS.se)), ])
} # }}}
