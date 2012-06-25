## get rid of common SNPs, chrX, chrY loci from a SummarizedExperiment 
dropUnclusterables <- function(SE, SNPs=NULL,RPTs=NULL,assay=NULL,byMad=FALSE){ 

  message('Excluding features on sex chromosomes...')
  unclusterable <- which(seqnames(rowData(SE)) %in% c('chrX','chrY'))

  if(unique(genome(rowData(SE))) == 'hg19' && is.null(SNPs)) {
    message("Automatically masking common SNPs from hg19...")
    require(FDb.FDb.UCSC.snp135common.hg19)
    SNPs <- features(FDb.FDb.UCSC.snp135common.hg19)
  }
  if(!is.null(SNPs)) {
    snp <- which(rownames(SE) %in% names(subsetByOverlaps(rowData(SE), SNPs)))
    unclusterable <- union(unclusterable, snp)
  } 
  if(unique(genome(rowData(SE))) == 'hg19' && is.null(RPTs)) {
    message("Automatically masking repeat regions from hg19...")
    require(BSgenome.Hsapiens.UCSC.hg19)
    ## this could stand to be in an FDb...
    RPTs <- do.call(c, lapply(seqlevels(Hsapiens), function(ch) {
      GRanges(ch, union(masks(Hsapiens[[ch]])$RM, masks(Hsapiens[[ch]])$TRF))
    }))
  }
  if(!is.null(RPTs)) {
    rpt <- which(rownames(SE) %in% names(subsetByOverlaps(rowData(SE), RPTs)))
    unclusterable <- union(unclusterable, rpt)
  } 
  if(byMad==TRUE) {
    require(matrixStats) 
    message('Excluding features with MAD of 0 or NaN...')
    if(is.null(assay)) assay = names(assays(SE, withDimnames=F))[[1]]
    values(rowData(SE))$mad <- rowMads(assays(SE, withDimnames=F)[[assay]])
    lowmad <- which(is.na(values(rowData(SE))$mad))
    unclusterable <- union(unclusterable, lowmad)
  }
  clusterable <- setdiff(seq_along(names(rowData(SE))), unclusterable)
  return(SE[ clusterable, ])

} 
