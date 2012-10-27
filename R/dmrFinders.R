dmrSE <- function(SE, groups, ...) { # {{{
  require('charm')
  message('This assumes your SE has already been properly annotated...')
  dmrFinder(eset=NULL, groups, p=assays(SE, withDimnames=F)$betas, 
            chr=as.character(seqnames(SE)), pos=start(SE), pns=rownames(SE),
            package='IlluminaHumanMethylation', ...)
} # }}}

dmrMatrix <- function(x, groups, what='betas', ...) { # {{{
  require('charm')
  if(!all(substr(rownames(x), 1, 2) == 'cg')) {
    stop('Your matrix of values must contain (ONLY) valid cgXXXX probe names')
  }
  require(FDb.InfiniumMethylation.hg19)
  browser()
  probes <- features(FDb.InfiniumMethylation.hg19)[rownames(x)]
  if(what=='betas') {
    dmrFinder(eset=NULL, groups, p=x, 
              chr=as.character(seqnames(probes)), 
              pos=start(probes), pns=names(probes),
              package='IlluminaHumanMethylation', ...)
  } else { 
    dmrFinder(eset=NULL, groups, l=x, 
              chr=as.character(seqnames(probes)), 
              pos=start(probes), pns=names(probes),
              package='IlluminaHumanMethylation', ...)
  }
} # }}}
