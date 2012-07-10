## this has simplified quite a bit
cpgWeights <- function(SE, decay=1000, verbose=FALSE) { # {{{

  require(Matrix)
  if(!identical(SE, sort(SE))) stop("Please sort() your SummarizedExperiment!")
  if(length(unique(seqnames(rowData(SE))))>1) stop("One chromosome at a time!")
  wts <- Matrix(0, ncol=nrow(SE), nrow=nrow(SE))
  diag(wts) <- 1
  rownames(wts) <- colnames(wts) <- rownames(SE)
  ol <- findOverlaps(rowData(SE), resize(rowData(SE),2*decay,fix='center'))

  idx <- split(subjectHits(ol), queryHits(ol))

  ## task: figure out which adjacent elements are > decay bp apart.
  blocked <- as.numeric(which(elementLengths(idx)>1))

  ## there's probably a faster way to do this.
  for(i in blocked) { ## this should go into C or C++
    if(verbose) { # {{{ progress meter
      msg <- paste('Computing weights for feature', i, 'of', length(blocked))
      message(paste0(msg,' (', round(100*(i/length(blocked))), '%)...'))
    } # }}}
    nonzero <- idx[[i]]
    z <- start(rowData(SE)[idx[[i]]])
    subwts <- 1 - log(1 + abs(sapply(z, '-', z)), base=decay)
    wts[ nonzero, nonzero ] <- subwts
  }

  ## add block structure attribute
  attr(wts, 'blocked') <- blocked

  ## to get the loci: rownames(wts)[as.numeric(attr(wts,'blocked'))]
  return(wts)

} # }}}

## weight NAs as zeros and normalize; results in sparser, per-sample weights
adjWts <- function(vec, wts) { # {{{

  require(matrixStats)
  empty <- which(is.na(vec))
  blocked <- attr(wts, 'blocked')
  if(length(intersect(empty,blocked))>0) {
    for(i in intersect(empty, blocked)) {
      wts[i, ] <- wts[ ,i] <- 0
      wts[i,i] <- 1
    }
    ## renormalize the matrix after knocking out NA rows/columns
    wts <- Matrix(apply(wts, 2, function(x) return(x / sum(x)))) 
  }
  return(wts)

} # }}}

## FIXME: do the whole thing via block diagonal matrix multiplications
smoothSites <- function(tmp, wts, keep.NAs=FALSE, plotMe=FALSE) { # {{{
  
  emptyCells <- is.na(tmp)
  breakpoints <- setdiff(seq_len(nrow(wts)), blocked) ## enables blocking

  ## if requested, plot 'before' 
  if(plotMe==TRUE) plot(seq_len(nrow(tmp)), tmp[,1], lwd=3, lty=3, col=1)

  ## this is where we should do block-diagonal multiplication; farm out to C++?
  message('Block diagonal multiplication would greatly speed up the following:')
  for(i in 1:ncol(tmp)) {

    ## 'before' plot
    if(plotMe==TRUE) lines(seq_len(nrow(tmp)), tmp[,i], lwd=3, lty=3, col=i)

    ## exponentially smooth the assay values in column 'i':
    tmp[,i] <- as.numeric(tmp[,i] %*% adjWts(tmp[,i], wts))

    ## 'after' plot
    if(plotMe==TRUE) lines(seq_len(nrow(tmp)), tmp[,i], lwd=3, lty=1, col=i)

  }

  if(keep.NAs) is.na(tmp) <- emptyCells
  return(tmp)

} # }}}

## sparse matrix multiplication for exponential smoothing
cpgSmooth <- function (SE, assay=NULL, wts=NULL, decay=1000, keep.NAs=F){ # {{{

  GR <- rowData(SE)
  if(!identical(GR, sort(GR))) stop('Please sort() your SummarizedExperiment!')
  if(length(unique(seqnames(rowData(SE)))) > 1) stop('Please split by chrom!')
  if(is.null(wts)) wts <- cpgWeights(SE, decay=decay)
  wts <- Matrix(apply(wts, 2, function(x) return(x / sum(x)))) ## normalize
  print(paste('Smoothing', unique(seqnames(rowData(SE))), 
              paste0('(', nrow(SE), 'features)...')))
  smoothed <- smoothSites(asy.fast(SE, assay), wts, keep.NAs=keep.NAs)
  rownames(smoothed) <- rownames(SE)
  colnames(smoothed) <- colnames(SE)
  assays(SE)[[ ifelse(is.null(assay), 1, assay) ]] <- smoothed
  return(SE)

} # }}}
