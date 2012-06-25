sdSdMax <- function(SE, assay=NULL) {
  require(matrixStats)
  if(is.null(assay)) assay = names(assays(SE, F))[[1]]
  x <- assays(SE, F)[[assay]]
  return(rowSds(x,na.rm=T)/sqrt(rowMeans(x,na.rm=T)/(1-(rowMeans(x,na.rm=T)))))
}

