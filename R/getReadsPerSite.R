## fetch reads, eventually normalized, from a BED or BAM file
getReadsPerSite <- function(SEorGR, filename, normalize=FALSE, totalReads=NULL){
  if(is(SEorGR, 'SummarizedExperiment')) {
    GR <- rowData(SEorGR)
  } else {
    GR <- SEorGR
  }
  if(all(grepl('(bed|gz)$', filename, ignore=T))) {
    ## it's a bed file
    require(rtracklayer)
    message('You DID run bedextract reads.bed sitesOfInterest.bed, yes?')
    bed.GR <- import(filename, asRangedData=FALSE)
    values(GR)$reads <- countOverlaps(GR, bed.GR)
  } else if(grepl('bam$', filename, ignore=T)) {
    ## it's a BAM file
    stop('BAM support remains unfinished at the moment; try converting to BED.')
  }
  if(normalize==TRUE) {
    if(!is.null(totalReads)) {
      values(GR)$normalized <- values(GR)$reads/totalReads
    } else {
      values(GR)$normalized <- values(GR)$reads/length(bed.GR)
    }
    values(GR)$normalized <- values(GR)$normalized/sum(values(GR)$normalized)
  }
  message(paste('Overlapping read', ifelse(normalize, 'densities:', 'counts:'),
                'across', length(GR), 'elements:'))
  print(summary(values(GR)$reads))
  return(GR)
}
