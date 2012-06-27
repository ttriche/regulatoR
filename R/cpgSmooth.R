## this has simplified quite a bit
cpgWeights <- function(SE, decay=1000) { # {{{
  require(Matrix)
  if(length(unique(seqnames(rowData(SE))))>1) stop("One chromosome at a time!")
  wts <- Matrix(0, ncol=nrow(SE), nrow=nrow(SE))
  diag(wts) <- 1
  rownames(wts) <- colnames(wts) <- rownames(SE)
  ol <- findOverlaps(rowData(SE), resize(rowData(SE), 2*decay, fix='center'))
  idx <- split(subjectHits(ol), queryHits(ol))
  for(i in seq_len(length(idx))) { ## this should go in C
    if(length(idx[[i]]) > 1) {
      msg <- paste('Computing weights for feature', i, 'of', length(idx))
      message(paste0(msg, ' (', round(100*(i/length(idx))), '%)...'))
      nonzero <- idx[[i]]
      z <- start(rowData(SE)[idx[[i]]])
      subwts <- 1 - log(1+abs(sapply(z, '-', z)), base=decay)
      wts[ nonzero, nonzero ] <- subwts
    }
  }
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
smoothSites <- function(tmp, wts, keep.NAs=FALSE) { # {{{
  emptyCells <- is.na(tmp)
  for(i in 1:ncol(tmp)) tmp[,i] <- as.numeric(tmp[,i] %*% adjWts(tmp[,i], wts))
  if(keep.NAs) is.na(tmp) <- emptyCells
  return(tmp)
} # }}}

## sparse matrix multiplication for exponential smoothing
cpgSmooth <- function (SE, assay=NULL, wts=NULL, decay=1000, keep.NAs=F){ #{{{
  assay <- chooseAssay(SE, assay)
  SE <- sort(SE) 
  GR <- rowData(SE)
  if(length(unique(seqnames(rowData(SE)))) > 1) stop('Please split by chrom!')
  if(is.null(wts)) wts <- cpgWeights(SE, decay=decay)
  print(paste('Smoothing', unique(seqnames(rowData(SE))), 
              paste0('(', nrow(SE), 'features)...')))
  smoothed <- smoothSites(asy.fast(SE, assay), wts, keep.NAs=keep.NAs)
  rownames(smoothed) <- rownames(SE)
  colnames(smoothed) <- colnames(SE)
  assays(SE)[[1]] <- smoothed
  return(SE)
} # }}}
