# distance-/eocg-weighted smoothing of 450k probes (correlation wise)
# load("~/Dropbox/SEs/CD34.METH.se.rda")
require(dlmodeler)

# first, we will want a representation of the correlation matrix across normals
# this could be the blood normals, all of TCGA normal tissues, or whatever...

# we want to take an X-base-pair rolling weighted average within r2 > 0.8 or so
# so rig up an Rle of the probes on each chromosomal arm and use runmeans...

# compare (ANOVA) using the 450k CD34/CD19/NEUT normals vs. Andrew Smith's BSseq


# now compare using the repeat percentage, SNP content, and so forth,
# with and without smoothing.  Can we soft-threshold these "bad" probes?


# How does the variance of each probe across samples compare to that of BS-seq?

# How does each probe across samples compare for Kalman vs. 100bp?
