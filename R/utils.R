## biweight midcorrelation for robust relationship detection
bicor.test <- function(x, y) { # {{{
  b = bicor(x,y)
  n = min(length(na.omit(x)), length(na.omit(y)))
  stat = abs(b * sqrt( (n-2)/(1-(b^2))))
  p.value = 1 - pt( stat, df=n-2 )
  pr = cor.test(x, y)
  r = pr$estimate
  res = c(bicor=b, b.t=stat, 
          r, pr$statistic, 
          b.p.value=p.value, r.p.value=pr$p.value)
  return(res)
} # }}}

## avoid endless name checking 
chooseAssay <- function(SE, assay=NULL) { # {{{
  return(ifelse(is.null(assay), names(assays(SE, withDimnames=F))[[1]], assay))
} # }}}

## avoid dimNames overhead by default
asy.fast <- function(SE, assay=NULL, withDimnames=FALSE) { # {{{
  assays(SE, withDimnames=withDimnames)[[ chooseAssay(SE, assay) ]]
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

## mostly for df2GR() and SNP finding
fromChr <- function(seqs, prefix='chr') { # {{{
  for(i in rev(seq_len(nchar(prefix)))) {
    seqs <- gsub(paste0('^', substr(prefix, 1, i)), '', seqs)
  }
  return(seqs)
} # }}}

## mostly for df2GR() and SNP finding
toChr <- function(seqs, prefix='chr') { # {{{
  paste0(prefix, fromChr(seqs, prefix))
} # }}}

## mostly for dropUnclusterables
byNAs <- function(SE, max=0.5, assay=NULL) { # {{{
  rowNAs <- rowSums(is.na(asy.fast(SE, assay)), na.rm=T)/ncol(SE)
  split(SE, ifelse(rowNAs > max, 'tooManyNAs', 'OK'))
} # }}}

## mostly for pyroPlot
byList <- function(SE, GRL, segmentBy='state', type='within') { # {{{
  if(!is(GRL, 'GRangesList')) {
    if(!(segmentBy %in% names(values(GRL)))) {
      stop(paste0('Please provide a value $', segmentBy, ' to segment on'))
    } else { 
      GRL <- split(GRL, as.vector(values(GRL)[[segmentBy]]))
    }
  }
  ol <- findOverlaps(rowData(SE), GRL, type=type) # default may dump features
  ols <- split(queryHits(ol), subjectHits(ol))
  names(ols) <- names(GRL)[ as.numeric(names(ols)) ]
  lapply(ols, function(rows) SE[ rows, ])
} # }}}

## for testing combine(SE, SE)
GEO.colData <- function(x) { # {{{
  colData(x)[ , grep('(characteristics|title)', names(colData(x)))]
} # }}}

## for testing combine(SE, SE)
min.colData <- function(x) { # {{{
  colData(x)[ , grep('(title|histology|age|celltype|status|gender)', 
                     names(colData(x))) ]
} # }}}
