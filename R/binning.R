require(Repitools) ## for genomeBlocks()

## genomeBlocks carves things up nicely:
## sl <- na.omit(seqlengths(rowData(SE)))
## blox <- genomeBlocks(sl, names(sl), blockSize)
##
## then apply by chromosome: 
##
## lapply(byChr(blox), function(x) findOverlaps(x, rowData(SE)))
## 
## and take trimmed means as needed (will need to move this to C)
##
##

## FIXME: add focal binning (eg. 100bp bins outward from the TSS or whatever)
binByFeature <- function(SE, features, assay=NULL, size=100, n=1, trim=0.01) {
  
  x <- asy.fast(SE, assay)
 

  featBinned <- lapply(split(features, values(features)$name),
                       featureBins, size=size)
  names(featBinned) <- c(names(genes), names(mirs))
  browser()
  summarizeBinned(SE, featBinned, assay=assay)    
}

## equal sized binned blocks for regression modeling e.g. exprs ~ meth + CNV
featureBins <- function(x,size=100,up=1500,down=100,body=F,enh=T,dist=250000){
  y <- resize(x, size, fix='start') ## first bin 
  for(i in seq_len(up/size)) y <- c(flank(y[1], size), y)
  for(i in seq_len(down/size)) y <- c(y, flank(y[1], size, FALSE))
  if(body==TRUE) {
    while(end(range(y)) < end(range(x))) {
      y <- c(y,flank(y[length(y)],size,F))
    }
  }
  if(enh==TRUE) {   ## should make it so this can be supplied instead
    if(unique(na.omit(genome(x)))!='hg19') stop("Don't have enhancers for hg18")
    if(!exists('REFSEQ.TES.HG19')) data(REFSEQ.TES.HG19)
    up.enhs <- resize(subsetByOverlaps(REFSEQ.TES.HG19, flank(x,dist)), 
                      fix='center', size)
    dn.enhs <- resize(subsetByOverlaps(REFSEQ.TES.HG19, flank(x,dist,F)),
                      fix='center', size)
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

