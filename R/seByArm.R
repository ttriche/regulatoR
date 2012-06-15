seByArm <-
function(SE, GRorGRL=NULL, build='hg19') {
  stopifnot(is(SE, 'SummarizedExperiment'))
  if(!is.null(GRorGRL)) {
    stopifnot(unique(genome(GRorGRL)) == build)
    stopifnot(class(GRorGRL) %in% c('GRanges','GRangesList'))
    stopifnot(unique(genome(rowData(SE))) == unique(genome(GRorGRL)))
  } else {
    data(hg19.by.arm)
    GRorGRL <- hg19.by.arm
  }
  if(is(GRorGRL, 'GRanges')) GRL <- split(GRorGRL) else GRL <- GRorGRL
  lapply(GRL, function(x) SE[ names(subsetByOverlaps(rowData(SE), x)), ])
}
