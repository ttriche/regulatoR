## so that I can send this around:
df2GR <- function(df, keepColumns=FALSE, ignoreStrand=FALSE) { ## {{{
  stopifnot(class(df) == "data.frame")
  subs <- list(chromStart='start', chromEnd='end', chrom='chr', seqnames='chr')
  for(s in names(subs)) names(df) = gsub(s, subs[[s]], names(df), ignore=TRUE)  
  stopifnot(all(c("start", "end") %in% names(df)))
  if('genome' %in% names(attributes(df))) g <- attr(df, 'genome') else g <- NULL
  if(substr(df$chr, 1, 3)[1] != 'chr') df$chr <- paste('chr', df$chr, sep='')
  df <- subset(df, !is.na(start) & !is.na(end))
  if(!ignoreStrand && ("strand" %in% names(df))) {
    if(is.numeric(df$strand)) df$strand <- strandMe(df$strand)
    GR <- with(df, GRanges(chr, IRanges(start=start, end=end), strand=strand))
  } else {
    GR <- with(df, GRanges(chr, IRanges(start=start, end=end)))
  }
  if('name' %in% names(df)) {
    names(GR) <- df$name
    df$name <- NULL
  } else {
    names(GR) <- rownames(df)
  }
  if(keepColumns) {
    skipped = c("rangename","chr","start","end","width","strand")
    elementMetadata(GR) <- as(df[, setdiff(names(df), skipped), drop=F], 
                              "DataFrame")
  }
  if('X' %in% names(elementMetadata(GR))) {
    if(all(is.na(GR$X))) {
      GR$X <- NULL
    } else {
      names(elementMetadata(GR))[which(names(elementMetadata(GR))=='X')]='score'
    }
  }
  if(!is.null(g)) genome(GR) <- g
  return(GR)
} # }}}

require(GenomicRanges)
setAs("MIAME", "SimpleList",
  function(from) { # {{{
    to = list()
    for(i in slotNames(from)) if(i != '.__classVersion__') to[[i]]=slot(from, i)
    return(SimpleList(to))
  }
) # }}}

## this will eventually require a better approach
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
} # }}}

setAs("ExpressionSet", "SummarizedExperiment", function(from) eSetToSE(from)) 

## for MethyLumiSets (sorta awkward as GEO 450k data has different assay names)
if(require(methylumi)) { ## this is going into MethyLumi anyhow
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
