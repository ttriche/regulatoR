## this has simplified quite a bit by going block-wise
cpgWeights <- function(SE, decay=1000, verbose=FALSE) { # {{{

  require(Matrix)
  if(!identical(SE, sort(SE))) stop("Please sort() your SummarizedExperiment!")
  if(length(unique(seqnames(rowData(SE))))>1) stop("One chromosome at a time!")
  wts <- Matrix(0, ncol=nrow(SE), nrow=nrow(SE))
  rownames(wts) <- colnames(wts) <- rownames(SE)
  diag(wts) <- 1
  blocks <- IRanges(
    start = c(1, which(diff(start(rowData(SE))) >= decay) + 1),
    end = c(which(diff(start(rowData(SE))) >= decay), length(rowData(SE)))
  )
  blocks = blocks[ width(blocks) > 1 ] 

  ## FIXME: profile and optimize this (use foreach here instead?!)
  for(i in seq_len(length(blocks))) { ## this should go into C or C++
    if(verbose) { # {{{ ghetto progress meter
      msg <- paste('Computing weights for block', i)
      message(paste0(msg, ' (', floor(100*(i/length(blocks))), '% done)...'))
    } # }}}
    nonzero <- c(start(blocks[i]):end(blocks[i]))
    z <- start(rowData(SE)[nonzero])
    subwts <- 1 - log(1 + abs(sapply(z, '-', z)), base=decay)
    wts[ nonzero, nonzero ] <- subwts
  }

  ## retain block boundaries...
  attr(wts, 'blocks') <- blocks
  return(wts)

} # }}}

## weight NAs as zeros and normalize; results in sparser, per-sample weights
adjWts <- function(vec, wts) { # {{{

  empty <- which(is.na(vec))
  blocks <- attr(wts, 'blocks')
  if(length(empty) > 0) { ## FIXME: farm this out to C/C++
    holes <- IRanges(start=empty, width=rep(1, length(empty)))
    toFix <- findOverlaps(blocks, holes)
    byBlock <- split(toFix, queryHits(toFix))
    for(i in seq_len(length(byBlock))) {
      for(j in subjectHits(byBlock[i])) {
        wts[j, ] <- wts[ ,j] <- 0
        wts[j,j] <- 1
      }
      submat <- seq( start(blocks[queryHits(byBlock[[i]])]), 
                     end(blocks[queryHits(byBlock[[i]])]) )
      tmp <- wts[submat, submat] ## local copy 
      ## renormalize the block after knocking out the NA rows/columns
      wts[submat, submat] <- matrix(apply(tmp, 2, function(x) return(x/sum(x))))
    }
  }
  return(wts)

} # }}}

## FIXME: do the whole thing via block diagonal matrix multiplications
smoothSites <- function(tmp, wts, keep.NAs=FALSE, plotMe=FALSE) { # {{{
  
  emptyCells <- is.na(tmp)
  blocks <- attr(wts, 'blocks')

  ## if requested, plot 'before' 
  if(plotMe==TRUE) plot(seq_len(nrow(tmp)), tmp[,1], lwd=3, lty=3, col=1)

  ## this is where we should do block-diagonal multiplication; farm out to C++?
  message('Block diagonal multiplication would greatly speed up the following:')

  ## FIXME: operate on blocks, as in adjWts 
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
