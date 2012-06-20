## for both ExpressionSet and Methy[Lumi][M][Set]s
setAs("MIAME", "SimpleList",
  function(from) { # {{{
    to = list()
    for(i in slotNames(from)) if(i != '.__classVersion__') to[[i]]=slot(from, i)
    return(SimpleList(to))
  }
) # }}}

## this will eventually require a better approach
setAs("ExpressionSet", "SummarizedExperiment",
  function(from) { # {{{
    if(grepl('^GPL', annotation(from))) {
      message('This coercion expects to find at least an Entrez ID via fData()')
      if( all( c('CHR','MAPINFO') %in% fvarLabels(from)) ) {
        fdat <- fData(from)[, c('ID','CHR','MAPINFO','MAPINFO') ]
        names(fdat) <- c('name','chrom','chromStart','chromEnd')
        fdat[, 3:4] <- apply(fdat[, 3:4], 2, as.numeric)
        row.dat <- df2GR(fdat)
      } else {
        stop(paste("Don't know how to handle platform", annotation(from)))
      }
    } else {
      chip <- annotation(from)
      require(paste(chip, "db", sep = "."), character.only=TRUE)
      starts <- unlist(lapply(mget(featureNames(from), 
                                   envir=get(paste0(chip, 'CHRLOC'))), 
                              function(x) x[[1]]))
      ends   <- unlist(lapply(mget(featureNames(from), 
                                   envir=get(paste0(chip, 'CHRLOCEND'))), 
                              function(x) x[[1]]))
      chrs   <- unlist(lapply(mget(featureNames(from), 
                                   envir=get(paste0(chip, 'CHR'))), 
                              function(x) x[[1]]))
      strand <- ifelse(sign(starts) == '-1', '-', '+')
      toGR <- data.frame(chrom=chrs, chromStart=abs(starts), chromEnd=abs(ends),
                         strand=strand, name=featureNames(from))
      flipped <- which(abs(toGR$chromStart)>abs(toGR$chromEnd)) 
      toGR[ flipped, c('chromEnd','chromStart') ] <- 
        abs(toGR[ flipped, c('chromStart','chromEnd') ])
      row.dat <- df2GR(toGR)
    }
    
    asy.dat <- SimpleList()
    asy.dat$exprs = assayDataElement(from, 'exprs')[names(row.dat), ]
    SummarizedExperiment(assays=asy.dat,
                         rowData=row.dat,
                         colData=as(pData(from), 'DataFrame'),
                         exptData=as(experimentData(from), 'SimpleList'))

  }) # }}}
      
if(require(methylumi)) { ## this is going into MethyLumi anyhow
  msetToSe <- function(from) { # {{{
    require(FDb.InfiniumMethylation.hg19) 
    chip = gsub('^IlluminaHumanMethylation','HM',gsub('k$','',annotation(from)))
    row.dat <- getPlatform(chip)
    asy.dat <- SimpleList()
    if(is(from, 'MethyLumiM')) {
      asy.dat$mvals = assayDataElement(from, 'exprs')[names(row.dat), ]
    } else if(is(from, 'MethyLumiSet')) {
      asy.dat$betas = assayDataElement(from, 'betas')[names(row.dat), ]
    }
    if( all( c('methylated','unmethylated') %in% assayDataElementNames(from)) ){
      asy.dat$total = assayDataElement(from, 'methylated')[names(row.dat), ] +
                      assayDataElement(from, 'unmethylated')[names(row.dat), ]
    }
    SummarizedExperiment(assays=asy.dat,
                         rowData=row.dat,
                         colData=as(pData(from), 'DataFrame'),
                         exptData=as(experimentData(from), 'SimpleList'))
  } # }}}
  setAs("MethyLumiSet", "SummarizedExperiment", function(from) msetToSE(from))
  setAs("MethyLumiM", "SummarizedExperiment", function(from) msetToSE(from))
} 

setAs("SummarizedExperiment", "GRanges", # just the first assay element for Gviz
  function(from) { # {{{
    message('SummarizedExperiment to GRanges retains ONLY the first assay...')
    GR = rowData(from)
    values(GR) = assays(from)[[1]]
    return(GR)
}) # }}}
