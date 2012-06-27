## I can't live without these, I'll go insane
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
