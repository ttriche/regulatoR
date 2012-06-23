## get rid of common SNPs, chrX, chrY loci from a SummarizedExperiment 
dropUnclusterables <- function(SE, SNPs=NULL,assay=NULL,RPTs=NULL,byMad=FALSE){ 

  if(is.null(assay)) assay = names(assays(SE, withDimnames=F))[[1]]
  unclusterable <- which(seqnames(rowData(SE)) %in% c('chrX','chrY'))

  if(byMad==TRUE) {
    require(matrixStats) 
    lowmad <- which(is.na(rowMads(assays(SE, withDimnames=F)[[assay]])))
    unclusterable <- union(unclusterable, lowmad)
  }
  if(!is.null(SNPs)) {
    snp <- which(rownames(SE) %in% names(subsetByOverlaps(rowData(SE), SNPs)))
    unclusterable <- union(unclusterable, snp)
  }
  if(!is.null(RPTs)) {
    rpt <- which(rownames(SE) %in% names(subsetByOverlaps(rowData(SE), RPTs)))
    unclusterable <- union(unclusterable, rpt)
  } 
  clusterable <- setdiff(seq_along(names(rowData(SE))), unclusterable)
  return(SE[ clusterable, ])

} 
