mergeAligned <- function(x, y, ...) { # {{{
  stopifnot(nrow(x$counts) == nrow(y$counts))
  stopifnot(identical(x$annotation$EntrezID, y$annotation$EntrezID))
  if (length(list(...)) > 0L) y = do.call(mergeAligned, list(y, ...))
  x$counts = cbind(x$counts, y$counts)
  x$targets = c(x$targets, y$targets)
  colnames(x$counts) = x$targets
  if('pvals' %in% names(x) && 'pvals' %in% names(y)) {
    x$pvals = cbind(x$pvals, y$pvals)
    colnames(x$pvals) = x$targets
  } else {
    warning('pvals not found for both sets; dropping...')
  }
  return(x)
} # }}}

## molds raw counts into a SummarizedExperiment and adds RPKM as an alternative
alignedToSE <- function(aligned, annotations=NULL, build='HG19') { # {{{

  require(Rsubread)
  asy.dat <- SimpleList()
  if('rpkm' %in% names(aligned)) asy.dat$rpkm <- aligned$rpkm
  else asy.dat$rpkm <- alignedToRPKM(aligned)
  asy.dat$counts <- aligned$counts
  if('pvals' %in% names(aligned)) asy.dat$pvals <- aligned$pvals
  
  if(!is.null(annotations)) { # user supplies their own annotations as with TCGA
    row.dat <- annotations
  } else if(build != 'HG19') {
    stop('You need to provide your own annotations for any genome besides hg19')
  } else if('ExonLength' %in% names(aligned$annotation)) {
    data('NCBI.EXONS.HG19')     # there's probably a better way to do this
    row.dat <- NCBI.EXONS.HG19
  } else if(nrow(aligned$annotation) > 22000) { # upper bound on hg19 lincRNAs
    data('NCBI.GENES.HG19')     # there's probably a better way to do this, too
    row.dat <- NCBI.GENES.HG19
    # omit genes w/ no annotated location
    NCBI.IDs = aligned$annotation$EntrezID 
    omit = which(!(NCBI.IDs %in% values(NCBI.GENES.HG19)$EntrezID))
    asy.dat <- endoapply(asy.dat, function(x) x[ -omit, ])
  } else {
    data('UCSC.LINCRNAS.HG19')    # there's definitely a better way to do this
    row.dat <- UCSC.LINCRNAS.HG19
  }

  coltmp <- DataFrame(ID=1:ncol(aligned$counts))
  x <- SummarizedExperiment(assays=asy.dat, rowData=row.dat, colData=coltmp)
  sampleNames(x) <- aligned$targets
  return(x)

} # }}}

## same thing but defaults to RPM instead of RPKM as the "alternate" measure
alignedToMiRnaSE <- function(aligned, annotations=NULL, build='HG19') { # {{{

  require(Rsubread)
  asy.dat <- SimpleList()
  if('rpm' %in% names(aligned)) {
    asy.dat$rpm <- aligned$rpm
  } else { 
    asy.dat$rpm <- alignedToRPM(aligned)
  } 
  asy.dat$counts <- aligned$counts
  
  if(!is.null(annotations)) {
    row.dat <- annotations
  } else if(build != 'HG19') {
    stop('You need to provide your own annotations for any genome besides hg19')
  } else {
    if(!exists('UCSC.MICRORNA.HG19')) data('UCSC.MICRORNA.HG19')
    row.dat <- UCSC.MICRORNA.HG19
  }

  coltmp <- DataFrame(ID=1:ncol(aligned$counts))
  x <- SummarizedExperiment(assays=asy.dat, rowData=row.dat, colData=coltmp)
  sampleNames(x) <- aligned$targets
  return(x)

} # }}}
