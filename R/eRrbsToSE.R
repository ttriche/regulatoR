eRrrbsToSE <- function(filename, genome='hg19') {
  require(GenomicRanges)
  if(!exists('df2GR')) source("df2GR.R")
  headers <- c('name','chr','start','strand','coverage','M','U')
  types <- list('character','character','integer','character',
                'integer','numeric','numeric') 
  names(types) <- headers
  foo <- as.data.frame(scan(filename, what=types, skip=1))
  foo$name = paste(foo$name, foo$strand, sep='.')
  foo$end = foo$start = as.integer(foo$start)
  foo$strand = as.factor(gsub('F','+', gsub('R','-', foo$strand)))
  foo$coverage = as.integer(foo$coverage)
  foo$M = round((as.numeric(foo$M) * foo$coverage)/100)
  foo$U = round((as.numeric(foo$U) * foo$coverage)/100)
  foo = foo[ , c('name','chr','start','end','strand','coverage','M','U') ]
  asys = SimpleList(methylated=as.matrix(round(foo$M/foo$U, 2), ncol=1), 
                    coverage=as.matrix(foo$coverage, ncol=1))
  rdat = df2GR(foo[,1:5])
  cdat = DataFrame(name=gsub('_myCpG.txt.gz','',filename), file=filename)
  SummarizedExperiment(assays=asys, rowData=rdat, colData=cdat)
}
