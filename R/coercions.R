require(GenomicRanges) ## for DataFrame, GenomicRanges, and so forth 

df2GR <- function(df, keepColumns=F, ignoreStrand=F, prefix='chr') { ## {{{

  if(class(df) == 'DataFrame') df <- as(df, 'data.frame')
  if(class(df) != 'data.frame') stop('df must be a data.frame or DataFrame')

  ## tidy up column names to coerce
  subs <- c(chr='chrom',
            seqnames='chrom',
            start='chromStart', 
            end='chromEnd')
  for(s in names(subs)) names(df) = gsub(s, subs[s], names(df), ignore=TRUE)
  if(!all(unique(subs) %in% names(df))) {
    stop('df must have columns chrom, chromStart, chromEnd to proceed')
  }

  ## assign genome pre-emptively if possible
  if('genome' %in% names(attributes(df))) {
    g <- attr(df, 'genome') 
  } else {
    g <- NULL
  }

  ## fix seqnames if necessary 
  if(!all(grepl(prefix, df$chrom))) {
    df$chrom <- toChr(df$chrom, prefix=prefix)
  }

  ## fuss about any missing data, which will be dropped presently 
  if( any(is.na(df$chromStart)) ) warning('Dropping ranges w/chromStart == NA')
  if( any(is.na(df$chromEnd)) ) warning('Dropping ranges w/chromEnd == NA')
  if( any(is.na(df$chrom)) ) warning('Dropping ranges w/chrom == NA')
  df <- subset(df, !is.na(chromStart) & !is.na(chromEnd) & !is.na(chrom))

  if(ignoreStrand == FALSE && ("strand" %in% names(df))) {
    if(is.numeric(df$strand)) {
      df$strand <- strandMe(df$strand)
    }
    GR <- with(df, GRanges(chrom, 
                           IRanges(start=chromStart, 
                                   end=chromEnd), 
                           strand=strand))
  } else {
    GR <- with(df, GRanges(chrom, 
                           IRanges(start=chromStart, 
                                   end=chromEnd)))
  }
  if('name' %in% names(df)) {
    names(GR) <- df$name
    df$name <- NULL
  } else {
    names(GR) <- rownames(df)
  }
  if(keepColumns) {
    skipped = c("rangename","chrom","chromStart","chromEnd","width","strand")
    elementMetadata(GR) <- as(df[, setdiff(names(df), skipped), drop=F], 
                              "DataFrame")
  }

  ## chintzy hack from the original methLab
  if('X' %in% names(elementMetadata(GR))) {
    if(all(is.na(GR$X))) {
      GR$X <- NULL
    } else {
      names(elementMetadata(GR))[which(names(elementMetadata(GR))=='X')]='score'
    }
  }

  ## assign genome to GR if known 
  if(!is.null(g)) genome(GR) <- g

  return(GR)
} # }}}
setAs("data.frame", "GenomicRanges", function(from) df2GR(from)) 
setAs("DataFrame", "GenomicRanges", function(from) df2GR(from)) 

setAs("MIAME", "SimpleList",
  function(from) { # {{{
    to = list()
    for(i in slotNames(from)) if(i != '.__classVersion__') to[[i]]=slot(from, i)
    return(SimpleList(to))
  }
) # }}}

## this will eventually require a better approach; add Illumina expr support 
eSetToSE <- function(from) { # {{{
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
    starts <- unlist(lapply(mget(featureNames(from), ifnotfound=NA,
                                 envir=get(paste0(chip, 'CHRLOC'))), 
                            function(x) x[[1]]))
    ends   <- unlist(lapply(mget(featureNames(from), ifnotfound=NA,
                                 envir=get(paste0(chip, 'CHRLOCEND'))), 
                            function(x) x[[1]]))
    chrs   <- unlist(lapply(mget(featureNames(from), ifnotfound=NA,
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
} # }}}
setAs("ExpressionSet", "SummarizedExperiment", function(from) eSetToSE(from)) 

## for MethyLumiSets (sorta awkward as GEO 450k data has different assay names)
if(require(methylumi)) { ## this is going into MethyLumi anyhow, eventually...
  msetToSE <- function(from) { # {{{
    require(FDb.InfiniumMethylation.hg19) 
    chip=gsub('^IlluminaHumanMethylation','HM',gsub('k$','',annotation(from)))
    row.dat <- getPlatform(chip)
    asy.dat <- SimpleList()
    features <- intersect(featureNames(from), names(row.dat))
    if(is(from, 'MethyLumiM')) {
      asy.dat$mvals = assayDataElement(from, 'exprs')[features, ]
    } else if(is(from, 'MethyLumiSet')) {
      asy.dat$betas = assayDataElement(from, 'betas')[features, ]
    }
    if(all( c('methylated','unmethylated') %in% assayDataElementNames(from))){
      asy.dat$total = assayDataElement(from, 'methylated')[features, ] +
                      assayDataElement(from, 'unmethylated')[features, ]
    }
    SummarizedExperiment(assays=asy.dat,
                         rowData=row.dat[features],
                         colData=as(pData(from), 'DataFrame'),
                         exptData=as(experimentData(from), 'SimpleList'))
  } # }}}
  setAs("MethyLumiSet", "SummarizedExperiment", function(from) msetToSE(from))
  setAs("MethyLumiM", "SummarizedExperiment", function(from) msetToSE(from))
} 

## for Gviz
setAs("SummarizedExperiment", "GRanges", # just the first assay element for Gviz
  function(from) { # {{{
    message('SummarizedExperiment to GRanges retains ONLY the first assay...')
    GR = rowData(from)
    values(GR) = assays(from)[[1]]
    return(GR)
}) # }}}
