## sparse correlation weighted regression for regulation
## this is kind of the heart of the entire package, in fact
## 
## example: y ~ (3 * mean(c(x2, x3, x4, x5, x6))) - 2*x15 + 5*x20
##          X is a 100xN matrix with only the above predictors, the rest noise
##          x2, x3, x4, x5, and x6 are highly multicollinear and within 100bp 
##          x15 is an enhancer locus not correlated w/others and far away
##          x20 is the exonic mean copy number, relative to normal.
## 
## goal:    fuse highly correlated predictors (x2:x6, x14:x15) eg. by bucketing
##          select sparse set of predictors eg. by using sparsenet to choose
##          call epigenetic regulation based on the robust MODEL probability
## 
## example: 
##
## notes:   the MODEL test statistic is used to gauge significance, 
##          so there must be some sort of a penalty in order to sparsify it,
##          but right now the best I can come up with is a two-step filter.
##
## details: just use the output from summarizeBinned() on CNV and 5MC 
##          to feed rlm( exprs ~ CNV + METH ) and get a model significance.
##
require(sparsenet)

