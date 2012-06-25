## FIXME: add focal binning (eg. 100bp bins outward from the TSS or whatever)
tssBin <- function(SE, bin=100, smoothme=F, assay=NULL, decay=1000, impute=T) {
  stopifnot(unique(genome(rowData(SE))) == 'hg19')
  if(is.null(assay)) assay <- names(assays(SE, FALSE))[[1]]
  if(smoothme) {
    x <- cpgSmooth(SE, assay=assay, decay=decay, impute=impute)
  } else {
    x <- assays(SE, F)[[assay]]
  }
  message('Need to add binning here: pull up TxDB and mean(findOverlaps(...))?')
  browser()
}
