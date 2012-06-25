cpgSmooth <- function (SE, assay = NULL, decay = 1000, impute = TRUE) {
  require(impute)
  require(parallel)
  options("mc.cores"=min(12, detectCores())) ## will thrash otherwise
  GR <- rowData(SE)[ grep('^cg', rownames(SE)), ] 
  if(is.null(assay)) assay = names(assays(SE, withDimnames = F))[[1]]
  stopifnot(unique(genome(rowData(SE))) == 'hg19')
  message('Subsetting by arm...')
  data(hg19.by.arm)
  GR.by.arm <- lapply(hg19.by.arm, function(x) subsetByOverlaps(GR, x))
  message(paste('Parallelizing by arm, over', options("mc.cores"), 'cores...'))
  by.arm <- mclapply(GR.by.arm, function(GR.arm) {
    print(paste('Locating neighbors on', unique(seqnames(GR.arm)), '...'))
    w <- findOverlaps(resize(GR.arm, 2 * decay, fix = "center"), GR.arm)
    idx <- split(w, queryHits(w))
    print(paste('Computing weights for', unique(seqnames(GR.arm)), '...'))
    getWeights <- function(x) { # {{{
      raw <- 1 - log(base = decay,
                     distance(GR.arm[queryHits(x)], GR.arm[subjectHits(x)]) + 1)
      return(raw / sum(raw))
    } # }}}
    wts <- lapply(idx, getWeights)
    arm.probes <- match(names(GR.arm), rownames(SE))
    tmp <- assays(SE, F)[[assay]][arm.probes, ]
    if(impute == TRUE && anyMissing(tmp)) tmp <- impute.knn(tmp)$data
    print(paste('Smoothing', unique(seqnames(GR.arm), '...')))
    reWeight <- function(x) { # {{{
      t(t(tmp[subjectHits(idx[[x]]), , drop = F]) %*% wts[[x]])
    } # }}}
    smoothed <- do.call(rbind, lapply(as.numeric(names(idx)), reWeight))
    rownames(smoothed) <- names(GR.arm)
    return(smoothed)
  })
  return(by.arm)
}
