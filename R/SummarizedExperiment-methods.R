
## I can't live without these, I'll go insane
setMethod("$", "SummarizedExperiment", function(x, name) { # {{{
  return(colData(x)[[name, exact=FALSE]])
}) # }}}
setMethod("$<-", "SummarizedExperiment", function(x, name, value) { # {{{
  colData(x)[[ name ]] <- value
  return(x)
}) # }}}
setMethod("seqinfo", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  seqinfo(rowData(x))
}) # }}}
setMethod("genome", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  genome(rowData(x))
}) # }}}

## there will be more in the future, this I guarantee.
