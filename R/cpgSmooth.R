cpgSmooth <- function (SE, assay = NULL, decay = 1000, impute = TRUE) {

  library(parallel)
  data(hg19.by.arm)
  GR <- rowData(SE)[ grep('^cg', rownames(SE)), ] 
  if(is.null(assay)) assay = names(assays(SE, withDimnames = F))[[1]]
  message('Subsetting by arm...')
  GR.by.arm <- lapply(hg19.by.arm, function(x) subsetByOverlaps(GR, x))
  message('Parallelizing by arm...')
  by.arm <- mclapply(GR.by.arm, function(GR.arm) {
    print('Locating neighbors...')
    w <- findOverlaps(resize(GR.arm, 2 * decay, fix = "center"), GR.arm)
    idx <- split(w, queryHits(w))
    print('Computing weights...')
    wts <- lapply(idx, function(x) { # {{{
      raw <- 1 - log(base = decay, 
                     distance(GR.arm[queryHits(x)],GR.arm[subjectHits(x)])+1)
      return(raw/sum(raw)) # normalize
    }) # }}}
    arm.probes <- match(names(GR.arm), rownames(SE))
    tmp <- assays(SE, F)[[assay]][arm.probes, ]
    if (impute == TRUE && anyMissing(assays(SE, F)[[assay]])) { # {{{
      message('Imputing NAs...')
      require(impute)
      tmp <- impute.knn(tmp)$data
    } # }}}
    print('Smoothing...')
    smoothed <- do.call(rbind, lapply(as.numeric(names(idx)), function(x) {
      t(t(tmp[subjectHits(idx[[x]]), , drop = F]) %*% wts[[x]])
    }))
    rownames(smoothed) <- names(GR.arm)
    return(smoothed)
  })
  smoothed <- do.call(rbind, by.arm)
  return(smoothed[ match(names(GR), rownames(smoothed)), ])

}
