## fetch reads, eventually normalized, from a BED or BAM file
getReadsPerSite <- function(filename, SEorGR, normalize=TRUE) {
  if(is(SEorGR, 'SummarizedExperiment')) SEorGR <- rowData(SEorGR)
  if(grepl('(bed|gz)$', filename, ignore=T)) {
    ## it's a bed file
    require(rtracklayer)
    message('You could use bedops for this, and it would be much faster...')
    bed.GR <- import(filename, asRangedData=FALSE)
    values(SEorGR)$counts <- countOverlaps(SEorGR, bed.GR)
  } else if(grepl('bam$', filename, ignore=T)) {
    ## it's a BAM file
    stop('BAM support remains unfinished at the moment; try converting to BED.')
  }
  if(normalize==TRUE) {
    values(SEorGR)$counts <- values(SEorGR)$counts/length(bed.GR)
  }
  return(SEorGR)
}
