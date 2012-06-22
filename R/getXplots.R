getXplots<-function(x, nProbes=500, col.fun='jet', clinical=NULL){

	require('GMD')
	require('fastcluster')
	require('GenomicRanges')
  if(!is(x, 'SummarizedExperiment')) {
    stop('You need to provide a SummarizedExperiment for this to work')
  }
  if( !('gender' %in% names(colData(x))) && !is.null(clinical) ){ 
    colData(x)[['gender']] <- clinical$gender
  }
  XInd <- sample( names(split(rowData(x), seqnames(rowData(x)))$chrX), nProbes)
  ## FIXME: assign gender based on mean methylation (or whatever the assay is)
	if(!is.null(colData(x)$gender)) {
    colSide <- ifelse(is.na(colData(x)$gender), 'white',
                      ifelse(colData(x)$gender=='M', 'lightblue', 'pink'))
  } else {
    colSide = rep('white', ncol(x))
  }
  name = as.character(match.call()["x"])
  asy <- names(assays(x, withDimnames=F))[[1]]
	tmp <- assays(x)[[asy]][XInd, ]
  jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                            "yellow", "#FF7F00", "red", "#7F0000"))
  bw <- colorRampPalette(c('white','gray','black'))
	capture.output({
    clusts <- suppressWarnings(
                heatmap.3(tmp, scale="none", trace="none", 
                          ColIndividualColors=colSide, color.FUN=get(col.fun), 
                          dendrogram='none', kc=2, labCol=NULL, 
                          labRow=1:nProbes, Colv=T, Rowv=TRUE,
                          main=paste('chrX clustering for', name)))$col.clusters
  })
  clusts <- as.factor(clusts)
  message('Assigned chrX cluster (assuming this is DNA methylation data):')
  return(clusts)
}


