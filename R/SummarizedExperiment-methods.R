setMethod("$", "SummarizedExperiment", function(x, name) { # {{{
  return(colData(x)[[name, exact=FALSE]])
}) # }}}
setMethod("$<-", "SummarizedExperiment", function(x, name, value) { # {{{
  colData(x)[[ name ]] <- value
  return(x)
}) # }}}
setMethod("sort", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  x[ names(sort(rowData(x))), ] 
}) # }}}
setMethod("genome", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  genome(rowData(x))
}) # }}}
setMethod("seqinfo", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  seqinfo(rowData(x))
}) # }}}
setMethod("combine", signature=signature(x="SummarizedExperiment", y="SummarizedExperiment"), function(x, y, ...) { # {{{
              if (class(x) != class(y)) {
                stop(paste("Error: objects must be the same class, but are ",
                           class(x), ", ", class(y), sep=""))
              }
              if( all(is.na(genome(rowData(x)))) || 
                  all(is.na(genome(rowData(y)))) ||
                  unique(na.omit(genome(rowData(x)))) !=
                  unique(na.omit(genome(rowData(y)))) ) {
                stop("Error: x and y have differing or unspecified genomes")
              }
              ## FIXME: allow for "packing out" missing features using NAs
              if( length(intersect(rownames(x), rownames(y))) < nrow(x) || 
                  length(intersect(rownames(x), rownames(y))) < nrow(y) ) {
                stop("Error: x and y have differing features, cannot combine")
              }
              browser()
              commonAsys <- intersect(names(x@assays), names(y@assays))
              names(commonAsys) <- commonAsys
              if(length(commonAsys) < 1) stop('Error: no assays in common')
              combineAssay <- function(assay, x, y) {
                cbind( assays(x, withDimnames=F)[[assay]],
                       assays(y[rownames(x), ], withDimnames=F)[[assay]] )
              }
              SummarizedExperiment(
                assays=lapply(commonAsys, combineAssay, x=x, y=y),
                colData=merge(colData(x), colData(y), all=TRUE),
                rowData=rowData(x)
              )
          }) # }}}
