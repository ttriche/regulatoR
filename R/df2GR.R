df2GR <- function(df, keepColumns=FALSE, ignoreStrand=FALSE) { # orig. from KDH 
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
}
