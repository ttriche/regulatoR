## FIXME: add focal binning (eg. 100bp bins outward from the TSS or whatever)
binByFeature <- function(SE, size=100, smooth=F,assay=NULL,decay=1000,impute=T){

  stopifnot(unique(genome(rowData(SE))) == 'hg19')
  if(is.null(assay)) assay <- names(assays(SE, FALSE))[[1]]
  if(smoothme) x <- cpgSmooth(SE, assay=assay, decay=decay, impute=impute)
  else x <- assays(SE, F)[[assay]]

  require(parallel)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)

  mirs <- microRNAs(TxDb.Hsapiens.UCSC.hg19.knownGene)
  names(mirs) <- values(mirs)$mirna_id
  values(mirs) <- DataFrame(name=names(mirs))
  mirs <- split(mirs, names(mirs))

  require(org.Hs.eg.db)
  genes <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 'gene')
  names(genes) <- mget(names(genes), org.Hs.egSYMBOL, ifnotfound=NA)
  
  features <- c(genes, mirs, .ignoreElementMetadata=TRUE) 
  values(features)$name <- c(lapply(names(genes), function(x) {
                               rep(x, length(genes[[x]]))
                             }), names(mirs))
  featBinned <- lapply(split(features, values(features)$name),
                       featureBins, size=size)
  names(featBinned) <- c(names(genes), names(mirs))
  browser()
  summarizeBinned(SE, featBinned, assay=assay)    
}

## equal sized binned blocks for regression modeling e.g. exprs ~ meth + CNV
featureBins <- function(x, size=100, up=1500, body=F, enhancers=T, dist=250000){
  y <- resize(x, size, fix='start') ## first bin 
  for(i in 1:(up/size)) y <- c(flank(y[[1]], size), y)
  if(body==TRUE) { 
    while(end(range(y)) < end(range(x))) y <- c(y,flank(y[[length(y)]],size,F))
  }
  if(enhancers==TRUE) {   ## should make it so this can be supplied instead
    data(REFSEQ.TES.HG19)
    up.enhs <- resize(subsetByOverlaps(REFSEQ.TES.HG19, flank(x,dist)), size)
    dn.enhs <- resize(subsetByOverlaps(REFSEQ.TES.HG19, flank(x,dist,F)), size)
    y <- c(up.enhs, y, dn.enhs)
  }
  return(y)
} 

## this function will feed directly into localRegression()
summarizeBinned <- function(SE, binning, assay=NULL) {
  if(is.null(assay)) assay <- names(assays(SE,F))[[1]]
  if(!is(binning, 'GRangesList')) binning <- split(binning)
  lapply(binning, function(x) { 
    sapply(x, function(y) {
      colMeans(assays(SE[which(rownames(SE) %in% 
                               names(subsetByOverlaps(rowData(SE), y))), ], 
                      withDimnames=FALSE)[[assay]], na.rm=TRUE)
    })
  })
}

