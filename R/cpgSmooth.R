## this has simplified quite a bit
cpgWeights <- function(SE, decay=1000) { # {{{
  require(Matrix)
  if(length(unique(seqnames(rowData(SE))))>1) stop("One chromosome at a time!")
  z = start(rowData(SE))
  wts = Matrix(1-log(sapply(z, function(x) pmin(abs(x-z)+1,decay)), base=decay))
  rownames(wts) <- rownames(SE)
  colnames(wts) <- rownames(SE)
  return(wts)
} # }}}

## weight NAs as zeros and normalize; results in sparser, per-sample weights
adjWts <- function(vec, wts) { # {{{
  require(matrixStats)
  empty <- which(is.na(vec))
  if(any(empty)) wts[empty, ] <- wts[, empty] <- 0
  Matrix(apply(wts, 2, function(x) {
    div = sum(x)
    if(div == 0) return(x)
    else return(x / div)
  }))
} # }}}

## FIXME: do the whole thing via matrices
smoothSites <- function(tmp, wts) {
  emptyCells <- is.na(tmp)
  for(i in 1:ncol(tmp)) tmp[,i] <- as.numeric(tmp[,i] %*% adjWts(tmp[,i], wts))
  is.na(tmp) <- emptyCells
  return(tmp)
}

## sparse matrix multiplication for exponential smoothing
cpgSmooth <- function (SE, assay=NULL, w=NULL, decay=1000) {
  assay <- chooseAssay(SE, assay)
  SE <- sort(SE) 
  GR <- rowData(SE)
  if(length(unique(seqnames(rowData(SE)))) > 1) stop('Please split by chrom!')
  if(is.null(w)) wts <- cpgWeights(SE, decay=decay)
  print(paste('Smoothing', unique(seqnames(rowData(SE))), 
              paste0('(', nrow(SE), 'features)...')))
  smoothed <- smoothSites(asy.fast(SE, assay), wts)
  rownames(smoothed) <- rownames(SE)
  colnames(smoothed) <- colnames(SE)
  assays(SE)[[1]] <- smoothed
  return(SE)
}
