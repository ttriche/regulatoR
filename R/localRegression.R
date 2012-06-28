## sparse correlation weighted regression for regulation
## this is kind of the heart of the entire package
##
## setup:   y ~ (3 * mean(c(x2, x3, x4, x5, x6))) - 2*x15 + 5*x20
##          X is a 100xN matrix with only the above predictors, the rest noise
##          x2, x3, x4, x5, and x6 are highly multicollinear and within 100bp 
##          x15 is an enhancer locus not correlated w/others and far away
##          x20 is the exonic mean copy number, relative to normal.
## 
## method:  fit a sparse LTS model y ~ x[ raw | binned | smooth | +cnv? ]
##          fit a robust lm on nonzero coefficients (y ~ x[ , coef(lts)>0 ])
##          fit a robust lm on the intercept (y ~ 1)
##          obtain a p-value for the improvement from including coefficients
##          this p-value is then adjusted genomewide to rank 'silencing' events
##
##          optional: use dp-means to visualize the silencing calls. 
## 
## example: fit = localRegression(SE, 'ECDH', 'binned', 100)
##
library(robust)
library(robustHD)

