## distance-/decay-weighted smoothing of 450k probes (correlation wise)
##
## this will need to be turned into a compiled function with 99.9% probability
##
distanceWeightedSmooth <- function(SE, assay=NULL, decay=1000) {

  library(parallel)
  GR <- rowData(SE)
  if(is.null(assay)) assay = names(assays(SE, withDimnames=F))[[1]]
  w <- findOverlaps(resize(GR, 2*decay, fix='center'), GR)
  idx <- split(w, queryHits(w))
  wts <- mclapply(idx, function(x) { # decaying correlation:
    raw <- 1 - log(base=decay, distance(GR[queryHits(x)], GR[subjectHits(x)])+1)
    return(raw/sum(raw)) # normalized weights
  }) 
  smoothed <- do.call(rbind, mclapply(as.numeric(names(idx)), function(x) {
    t(t(assays(SE, F)[[assay]][ subjectHits(idx[[x]]), , drop=F]) %*% wts[[x]])
  }))
  message('This function does not yet deal with NAs -- beware.  Or impute.')
  message('Seriously -- distance-weighted smoothing needs testing and tuning.') 
  return(smoothed)

}

# compare (ANOVA) using the 450k CD34/CD19/NEUT normals vs. Andrew Smith's BSseq
# also compare to the normal bone marrows from Melnick & Figueroa (eRRBS)
# and perhaps to PBMNCs from Esteller (since I'm going to use them as normals!)
# now compare using the repeat percentage, SNP content, and so forth,
# with and without smoothing.  Can we soft-threshold these "bad" probes?
