## avoid endless name checking 
chooseAssay <- function(SE, assay=NULL) { # {{{
  return(ifelse(is.null(assay), names(assays(SE, withDimnames=F))[[1]], assay))
} # }}}

## avoid dimNames overhead by default
asy.fast <- function(SE, assay=NULL) { # {{{
  assays(SE, withDimnames=FALSE)[[ chooseAssay(SE, assay) ]]
} # }}}

## mostly for CpG smoothing and silencing calls
byChr <- function(x) { # {{{
  if(is(x, 'SummarizedExperiment')) {
    split(x, as.vector(seqnames(rowData(x))))
  } else if(is(x, 'GenomicRanges')) {
    split(x, seqnames(x))
  } else if(is(x, 'List')) {
    split(x, x$chr, drop=T)
  }
} # }}}

## mostly for dropUnclusterables
byNAs <- function(SE, max=0.5, assay=NULL) { # {{{
  rowNAs <- rowSums(is.na(asy.fast(SE, assay)), na.rm=T)/ncol(SE)
  split(SE, ifelse(rowNAs > max, 'tooManyNAs', 'OK'))
}

## mostly for pyroPlot
byList <- function(SE, GR) {
  if(!'state' %in% names(values(GR))) stop('Need a $state value to segment on')
  else GRL <- split(GR, as.vector(values(GR)$state))  # should be a factor
  ol <- findOverlaps(rowData(SE), GRL, type='within') # may dump some probes
  
  message('Not finished...')
} 

## for testing cpg smoothing
fakeProbes <- function(x) {
  matrix(rbeta(x^2, 2, 6.7), ncol=x)
}

## as above
fakeWeights <- function(x) {
  require(Matrix)
  Matrix(cor(matrix(rbeta(x^2, 2, 6.7), ncol=x)))
}

## as above
fakeNAs <- function(tmp, prop=0.2) {
  is.na(tmp) <- matrix(rbinom(ncol(tmp)*nrow(tmp), 1, prop), ncol=ncol(tmp))==1
  tmp
}
