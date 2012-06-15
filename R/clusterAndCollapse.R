clusterAndCollapse <-
function(x, minSep, G=c(1:5)) {  ## this should recurse
  require(mclust)
  results = list()
  results$retain = 0
  results$values = x
  mcl = Mclust(x, G=G) # look for small i# of groups by default
  means = mcl$parameters$mean
  diff.mat = sapply(means, function(x) sapply(means, '-', x))
  if(all(diff.mat < minSep)){
    return(results)
  } else if(any(is.na(x))) {                                                    
    message('NA/NaN values found, skipping... impute if you want to use these')
    return(results)
  } else if(any(diff.mat[ -which(diff.mat==0) ] < minSep)) { 
    rounds = 1
    i = rev(order(table(as.factor(mcl$classification))))[1]
    while(rounds < mcl$G) { # {{{ the business: collapse if sep < minSep 
      if(length(which(mcl$classification == i)) == 0) break()
      newmean = means[i]
      groups = which(table(as.factor(mcl$classification)) > 0)
      means = tapply(x, mcl$classification, mean, na.rm=T)
      diffs = sapply(means, function(x) abs(x - newmean))
      names(diffs) = names(means)
      if(all(diffs > minSep)) break()
      rounds = rounds + 1 
      to.merge = groups[ which(diffs[groups] < minSep) ]
      mcl$classification[ which(mcl$classification %in% to.merge) ] = i
      means = tapply(x, mcl$classification, mean, na.rm=T)
      if(length(means) < 2) break()
      newmean = mean(x[ which(mcl$classification == i) ], na.rm=T)
      diffs = sapply(means, function(x) abs(x - newmean))
      if(all(diffs > minSep)) break()
      groupnames = names(rev(sort(diffs*table(mcl$classification))))
      for(j in seq_along(groupnames)) {
        if(groupnames[j] == i) break()
        i = j
        if(diffs[groupnames[i]] > minSep) break()
      } 
    } # }}}
    means = tapply(x, mcl$classification, mean, na.rm=T)
    results$values = means[ as.character(mcl$classification) ]
    results$retain = length(unique(results$values))
    return(results)
  }
}
