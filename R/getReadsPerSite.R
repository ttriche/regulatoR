## fetch reads, eventually normalized, from a BED or BAM file
getReadsPerSite <- function(filename, SEorGR) {
  if(is(SEorGR, 'SummarizedExperiment')) SEorGR <- rowData(SEorGR)
  if(grepl('(bed|gz)$', filename, ignore=T)) {
    ## it's a bed file
    require(rtracklayer)
    message('You could use bedops for this, and it would be much faster...')
    bed.GR <- import(filename)
    reads <- countOverlaps(SEorGR, bed.GR)
  } else if(grepl('bam$', filename, ignore=T)) {
    ## it's a BAM file
    message('BAM files can be real pigs, maybe convert it to a BED file first?')
    stop('BAM support remains unfinished at the moment')  
  }
  return(reads)
}
