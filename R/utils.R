## I use this all over now
chooseAssay <- function(SE, assay=NULL) {
  if(is.null(assay)) return(names(assays(SE, F))[[1]])
  else return(assay)
}

## mostly for CpG smoothing and silencing calls
byChr <- function(SE) split(SE, as.vector(seqnames(rowData(SE))))

## mostly for dropUnclusterables
byNAs <- function(SE, max=0.5, assay=NULL) {
  assay <- chooseAssay(SE, assay)
  rowNAs <- rowSums(is.na(assays(SE,F)[[assay]]))/ncol(SE)
  split(SE, ifelse(rowNAs > max, 'tooManyNAs', 'OK'))
}

## mostly for pyroPlot
byList <- function(SE, GR) {
  if(!'state' %in% names(values(GR))) stop('Need a $state value to segment on')
  else GRL <- split(GR, as.vector(values(GR)$state))  # should be a factor
  ol <- findOverlaps(rowData(SE), GRL, type='within') # may dump some probes
   
  message('Not finished...')
} 
