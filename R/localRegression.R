## sparse correlation weighted regression for regulation
## this is kind of the heart of the entire package
##
## setup:   y ~ (3 * mean(c(x2, x3, x4, x5, x6))) - 2*x15 + 5*x20
##          X is a 100xN matrix with only the above predictors, the rest noise
##          x2, x3, x4, x5, and x6 are highly multicollinear and within 100bp 
##          x15 is an enhancer locus not correlated w/others and far away
##          x20 is the exonic mean copy number, relative to normal.
## 
## method:  bin using FeatureBlocks, add 'TypeI' and 'TypeII' interactions, CNV?
##          select predictors using model <- FWDselect::qselection(crit='R2')
##          fit a robust lm on model$selection[ which.max(model$R2) ]
##          fit a robust lm on the intercept (y ~ 1)
##          obtain a p-value for the improvement from including coefficients
##          this p-value is then adjusted genomewide to rank 'silencing' events
##
## example: fit = localRegression(LAML, LAML.gene, 'ECDH', 'bin', 100)
##
localRegression <- function(mSE, eSE, gene, how='bin', binSize=100, inBins=F) {

  require(MASS)
  require(robust)
  require(FWDselect)

  if(!inBins) {

    ## create featureBlocks for the provided gene from the provided 5mC SE


  }

}
