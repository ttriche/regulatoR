require(Biobase)
require(GenomicRanges)

# SummarizedExperiment is a nice replacement for ExpressionSet and friends
setMethod("$", "SummarizedExperiment", function(x, name) { # {{{
  return(colData(x)[[name, exact=FALSE]])
}) # }}}
setMethod("$<-", "SummarizedExperiment", function(x, name, value) { # {{{
  colData(x)[[ name ]] <- value
  return(x)
}) # }}}
setMethod("pData", "SummarizedExperiment", function(object) { # {{{
  colData(object)
}) # }}}
setMethod("pData<-","SummarizedExperiment",function(object, value) { # {{{
  colData(object) <- value
  return(object)
}) # }}}
setMethod("fData", "SummarizedExperiment", function(object) { # {{{
  rowData(object)
}) # }}}
setMethod("fData<-","SummarizedExperiment",function(object, value) { # {{{
  rowData(object) <- value
  return(object)
}) # }}}
setMethod("sampleNames", "SummarizedExperiment", function(object) { # {{{
  colnames(object)
}) # }}}
setMethod("sampleNames<-","SummarizedExperiment",function(object, value) { # {{{
  colnames(object) <- value
  return(object)
}) # }}}
setMethod("featureNames", "SummarizedExperiment", function(object) { # {{{
  names(rowData(object))
}) # }}}
setMethod("featureNames<-","SummarizedExperiment",function(object, value) {# {{{
  names(rowData(object)) <- value
  return(object)
}) # }}}
setMethod("seqinfo", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  seqinfo(rowData(x))
}) # }}}
setMethod("genome", signature(x="SummarizedExperiment"), function(x) { # {{{ 
  genome(rowData(x))
}) # }}}
