## get rid of common SNPs, chrX, chrY loci from a SummarizedExperiment 
dropUnclusterables <- function(SE, SNPs=NULL,RPTs=NULL,assay=NULL) {

  message('Dropping features with > 50% NAs...')
  rowNAs <- function(SE) 
  unclusterable <- which(seqnames(rowData(SE)) %in% c('chrX','chrY'))
  message('Excluding features on sex chromosomes...')
  unclusterable <- which(seqnames(rowData(SE)) %in% c('chrX','chrY'))
  if(unique(genome(rowData(SE))) == 'hg19' && is.null(SNPs)) {
    message("Automatically masking common SNPs from hg19...")
    require(FDb.UCSC.snp135common.hg19)
    SNPs <- features(FDb.UCSC.snp135common.hg19)
    detach(package:FDb.UCSC.snp135common.hg19, unload=TRUE)
    gc(,T)
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
    detach(package:BSgenome.Hsapiens.UCSC.hg19, unload=TRUE) 
    gc(,T)
  }
  if(!is.null(RPTs)) {
    rpt <- which(rownames(SE) %in% names(subsetByOverlaps(rowData(SE), RPTs)))
    unclusterable <- union(unclusterable, rpt)
  } 
  clusterable <- setdiff(seq_along(names(rowData(SE))), unclusterable)
  return(SE[ clusterable, ])
} 
