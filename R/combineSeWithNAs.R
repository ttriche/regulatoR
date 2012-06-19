## Take a bunch of ragged SummarizedExperiment objects and inject NAs 
## for the ones which don't have observations at a feature where others do.
## This reduces to "find the superset of features and inject NAs elsewhere".
## The assumption is that some sort of smoothing will be done later on. 
## This was originally written to combine enhanced RRBS data for comparisons.
##
combineSeWithNAs <- function(x, y, ...) { 

  ## get some annoying preliminaries out of the way first...
  if(length(list(...)) > 0L) y <- do.call(combineSeWithNAs, list(y, ...))
  if(class(x)!=class(y)) stop(paste("class mismatch: x is a", class(x), ", ", 
                                     "while y is a", class(y), sep=""))
  if(names(values(x)) != names(values(y))) stop("names(values(rowData)) differ")
  else superset <- union(rowData(x), rowData(y))

  message("Adding NAs for unshared features...")

  stop('Ragged SE combination (eg. for eRRBS data) is not yet finished')

}
