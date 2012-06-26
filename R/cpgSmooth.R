getWeight <- function(x, GR, decay=1000) {
  raw = 1-log((distance(GR[queryHits(x)], GR[subjectHits(x)])+1), base=decay)
  names(raw) <- names(GR)[subjectHits(x)]
  raw / sum(raw) # normalized
}

## FIXME: do the whole thing via matrices
reWeightCpG <- function(x, tmp, idx, wts) {
  return(t(t(tmp[subjectHits(idx[[x]]), , drop = F]) %*% wts[[x]]))
}

## FIXME: use a large matrix, not a list
cpgWeight <- function (SE, decay=1000) {
  GR <- rowData(SE)
  ol <- findOverlaps(resize(GR, 2*decay, fix="center"), GR)
  idx <- split(ol, queryHits(ol))
  names(idx) <- names(GR)[ queryHits(ol) ]
  wts <- lapply(idx, getWeight, GR=GR, decay=decay)
  return(list(wts=wts, idx=idx))
}

## FIXME: use matrix multiplication to obtain the reweighted values
cpgSmooth <- function (SE, w=NULL, decay=1000, assay=NULL, impute=T) {
  require(impute)
  GR <- rowData(SE)
  if(is.null(assay)) assay = names(assays(SE, withDimnames = F))[[1]]
  if(length(unique(seqnames(rowData(SE)))) > 1) message('You CAN split by chr!')
  if(is.null(w)) w <- cpgWeight(SE, decay=decay)
  idx <- w$idx
  wts <- w$wts
  tmp <- assays(SE, F)[[assay]]
  if(impute == TRUE && anyMissing(tmp)) tmp <- impute.knn(tmp)$data
  print(paste('Smoothing', unique(seqnames(rowData(SE))), '...'))
  smoothed <- do.call(rbind, lapply(names(idx), 
                                    reWeightCpG, tmp=tmp, idx=idx, wts=wts))
  rownames(smoothed) <- rownames(SE)
  return(smoothed)
}
