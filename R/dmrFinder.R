
dmrFinder <- function(eset=NULL, groups, p=NULL, l=NULL, chr=NULL, pos=NULL, pns=NULL, sdBins=NULL, controlIndex=NULL, controlProbes=NULL, Indexes=NULL, filter=NULL, package=NULL, ws=7, verbose=TRUE, compare="all",  withinSampleNorm="loess", betweenSampleNorm="quantile", cutoff=0.995, sortBy="ttarea", removeIf=expression(nprobes<3), paired=FALSE, pairs=NULL, DD=NULL, COMPS=NULL, COMPS.names=NULL) { # {{{

  if(!is.null(removeIf) & !is.expression(removeIf)) stop("If not NULL, removeIf argument must be an expression.")
  groups = as.character(groups)
  if(paired & is.null(DD) & is.null(pairs)) stop("pairs argument must be provided if paired=TRUE.")
  if(paired & is.null(DD) & length(pairs)!=length(groups)) stop("pairs argument must be same length as groups.")

  if(!paired | is.null(DD)) { #then the compare arg will be used.
      if(identical(compare,"all")) compare=comp(groups)
      if(length(compare)%%2!=0) stop("compare must have an even number of elements.")
      if(length(cutoff)==1) cutoff <- rep(cutoff, length(compare)/2)
      if(length(compare)/2!=length(cutoff)) stop(length(compare)/2," comparisons requested but ", length(cutoff)," cutoff(s) provided.")
      if(!all(compare%in%groups)) stop("Not all groups specified in the compare argument are in the groups argument.")
  }

  args=list(filter=filter, ws=ws, betweenSampleNorm=betweenSampleNorm, 
	    withinSampleNorm=withinSampleNorm, sdBins=sdBins,
            cutoff=cutoff, sortBy=sortBy, compare=compare, 
            paired=paired, pairs=pairs, removeIf=removeIf)

  # dmrFinder must be given either eset or p/l,chr,pos,pns, and package.
  # If eset is supplied then all the latter will be taken from it (with any
  # that were given as arguments being ignored) (except p).
  # l=logit(p). Indexes=split(seq(along=pns),pns).
  if(is.null(eset)) {
      if (is.null("p") & is.null("l")) stop("p or l must be supplied.")
      Args = c("chr","pos","pns","package")  #"controlIndex" not needed
      nulls = sapply(Args,function(x) is.null(get(x)))
      if(any(nulls))
        stop(paste("The following arguments are missing:", paste(Args[nulls], collapse=", ")))
      lens = c( nrow(p), nrow(l), length(chr), length(pos), length(pns) ) 
      if(length(unique(lens))!=1)
        stop("p, l, chr, pos, and/or pns are incompatible.")
      stopifnot(length(groups)==max(ncol(p), ncol(l)))
	  chrpos <- paste(chr, pos)
	  dup <- duplicated(chrpos) | duplicated(chrpos, fromLast = TRUE)
	  index <- which(!is.na(chr) & !is.na(pos) & !is.na(pns) & !dup)
      index=index[order(chr[index],pos[index])]
      chr=chr[index]
      pos=pos[index]
      pns=pns[index]
#      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      if(!is.null(p)) p=p[index,]
      if(!is.null(l)) l=l[index,]
  } else if (is.character(eset)) {
          if (is.null("p") & is.null("l")) stop("p or l must be supplied.")
	  pdInfo=get(eset)
	  class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition, probeNames work
	  chr=pmChr(pdInfo)
	  pos=pmPosition(pdInfo)
	  chrpos <- paste(chr, pos)
	  dup <- duplicated(chrpos) | duplicated(chrpos, fromLast = TRUE)
	  if (!is.null(l)) {
		index=which(rowSums(is.na(l))==0 & !dup)	
        index=index[order(chr[index],pos[index])]
		l=l[index,]
	  } else {
		index=which(rowSums(is.na(p))==0 & !dup)	
        index=index[order(chr[index],pos[index])]
        p=p[index,]
	  }
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(pdInfo)[index]
#      if (is.null(controlIndex)) {
#          controlIndex <- which(getContainer(pdInfo) %in% controlProbes)
#       }
#      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      package = eset
      if(package=="pd.feinberg.mm8.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
      if(package=="pd.feinberg.hg18.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
  } else {
      stopifnot(length(groups)==ncol(eset))
      if(is.null(p) & is.null(l)){
		p <- methp(eset, 
					withinSampleNorm=withinSampleNorm, 
					betweenSampleNorm=betweenSampleNorm,
					controlIndex=controlIndex,
					controlProbes=controlProbes, 
					verbose=TRUE)
	  }
      chr=pmChr(eset)
      pos=pmPosition(eset)
	  chrpos <- paste(chr, pos)
	  dup <- duplicated(chrpos) | duplicated(chrpos, fromLast = TRUE)
	  if (!is.null(l)) {
		index=which(rowSums(is.na(l))==0 & !dup)	
        index=index[order(chr[index],pos[index])]
		l=l[index,]
	  } else {
		index=which(rowSums(is.na(p))==0 & !dup)			
        index=index[order(chr[index],pos[index])]
        p=p[index,]
	  }
      chr=chr[index]
      pos=pos[index]
      pns=probeNames(eset)[index]
#      if (is.null(controlIndex)) {
#          controlIndex=getControlIndex(eset, controlProbes=controlProbes)
#      }
#      controlIndex=which(index%in%controlIndex)
      if(!is.null(sdBins)) sdBins<-sdBins[index]
      package = annotation(eset)
      if(package=="pd.feinberg.mm8.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
      if(package=="pd.feinberg.hg18.me.hx1"){
      #add some code here to break up regions with gaps of >300 bp
      }
  }
  if(is.null(l)) {
	  l=log(p)-log(1-p)	
  }
  if(is.null(Indexes)) {
  	  Indexes=split(seq(along=pns),pns)
  }

  if(!paired) {
      tog = get.tog(l=l,groups=groups,compare=compare,verbose=verbose)
      lm=tog$lm
      ls=tog$ls
      ns=tog$ns
      COMPS=tog$COMPS
      tt = get.tt(lm=lm,ls=ls,ns=ns,COMPS=COMPS,Indexes=Indexes,filter=filter,ws=ws,verbose=verbose)
  } else {
      if(is.null(DD)) {
          DD0 = get.DD(l=l, groups=groups, compare=compare, pairs=pairs)
          DD    = DD0$DD
          COMPS = DD0$COMPS
          COMPS.names = DD0$COMPS.names
      } else { 
          if(is.null(COMPS)|is.null(COMPS.names)) stop("If DD is provided, COMPS and COMPS.names must also be provided.")
      }
      ttpaired = get.tt.paired(DD=DD,Indexes=Indexes,filter=filter,ws=ws,verbose=verbose)
      tt  = ttpaired$sT
      MD  = ttpaired$MD
      sMD = ttpaired$sMD
  }

  res=vector("list",ncol(tt))
  names(res)=colnames(tt)
  if(verbose) message("Finding DMRs for each pairwise comparison.")
  for(r in 1:nrow(COMPS)) {
      j = COMPS[r,1]
      k = COMPS[r,2]
      if(verbose) message("\n",colnames(tt)[r])
      if(!paired) DF=ifelse(ns[j]==1 & ns[k]==1, 1, ns[j]+ns[k]-2)
      if(paired)  DF=ifelse(ncol(DD[[r]])-1==0,1,ncol(DD[[r]])-1)
	
	  if (length(sdBins)==0) {
	      K=mad(tt[,r], na.rm=TRUE)*qt(cutoff[r],DF)	
	  }	else {
		  s <- tapply(tt[,r], sdBins, mad, na.rm=TRUE)
		  K=s[sdBins]*qt(cutoff[r],DF)	
	  }
      LAST=0
      segmentation=vector("numeric",nrow(tt))
      type=vector("numeric",nrow(tt))
      for(i in seq(along=Indexes)){
        if(verbose) if(i%%1000==0) message(".", appendLF=FALSE)
        Index=Indexes[[i]]
        y=tt[Index,r]
		if(length(sdBins)==0) {
			tmp=sign(y)*as.numeric(abs(y)>K)	
		} else {
			Ki <- K[Index]
			tmp=sign(y)*as.numeric(abs(y)>Ki)	
		}
        tmp2=cumsum(c(1,diff(tmp)!=0))+LAST
        segmentation[Index]=tmp2
        type[Index]=tmp
        LAST=max(tmp2)
      }
      
      Index = which(type!=0)
      rows  = length(unique(segmentation[Index]))
      res[[r]]=data.frame(
           chr=tapply(chr[Index],segmentation[Index],function(x) x[1]),
           start=tapply(pos[Index],segmentation[Index],min),
           end=tapply(pos[Index],segmentation[Index],max),
           p1=rep(NA, rows),
           p2=rep(NA, rows),
           regionName=tapply(pns[Index],segmentation[Index],function(x) x[1]),
           indexStart=tapply(Index,segmentation[Index],min),
           indexEnd=tapply(Index,segmentation[Index],max))
      length=res[[r]]$indexEnd-res[[r]]$indexStart+1
      res[[r]]$nprobes = length
## The following commented out part gives identical results to the uncommented out replacement code below if using l rather than p, and if odd number of samples in each group, when using p too, but otherwise lm's values are the means of the middle 2 values and transforming back from logit to p isn't exact so there will be very minor discrepancies.
#     if(!paired) {
#         if(is.null(p)) { #  We return log-ratios
#	      colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
#	      colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
#	      res[[r]]$m1=tapply(lm[Index,j],segmentation[Index],mean)
#             res[[r]]$m2=tapply(lm[Index,k],segmentation[Index],mean)
#             diff = res[[r]]$m1 - res[[r]]$m2
#             maxdiff=tapply(lm[Index,j]-lm[Index,k],segmentation[Index], 
#                            function(x) { x[which.max(abs(x))] })
#	  } else { # We return percentages
#             tmp1 = 1/(1+exp(-lm[Index,j]))
#             tmp2 = 1/(1+exp(-lm[Index,k]))
#	      res[[r]]$p1=tapply(tmp1,segmentation[Index],mean)
#             res[[r]]$p2=tapply(tmp2,segmentation[Index],mean)
#             diff = res[[r]]$p1 - res[[r]]$p2                
#	      maxdiff=tapply(tmp1-tmp2,segmentation[Index], 
#                            function(x) { x[which.max(abs(x))] })
#         }
#     } else {

      if(nrow(res[[r]])>0) {
          if(!paired) {
              grp1 = which(groups==colnames(lm)[j])
              grp2 = which(groups==colnames(lm)[k])
          } else {
              grp1 = which(groups==COMPS.names[r,1])
              grp2 = which(groups==COMPS.names[r,2])
          }
          if(is.null(p)) {
              colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
	      colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
              mat  = cbind(rowMedians(l[,grp1]), rowMedians(l[,grp2]))
              matr = t(apply(res[[r]][,c("indexStart","indexEnd")],1,
                             function(se) colMeans(mat[se[1]:se[2],,drop=FALSE])))
              res[[r]]$m1 = matr[,1]
              res[[r]]$m2 = matr[,2]
              diff = res[[r]]$m1 - res[[r]]$m2
          } else {
              mat  = cbind(rowMedians(p[,grp1]), rowMedians(p[,grp2]))
              matr = t(apply(res[[r]][,c("indexStart","indexEnd")],1,
                             function(se) colMeans(mat[se[1]:se[2],,drop=FALSE])))
              res[[r]]$p1 = matr[,1]
              res[[r]]$p2 = matr[,2]
              diff = res[[r]]$p1 - res[[r]]$p2
          }
          ## Max diff:
          matd = mat[,1]-mat[,2]
          maxdiff = apply(res[[r]][,c("indexStart","indexEnd")],1,
                          function(se) {
                              rt = matd[se[1]:se[2]]
                              rt[which.max(abs(rt))] })
     
          area   = abs(diff)*length
          ttarea = abs(tapply(tt[Index,r],segmentation[Index],mean)) *length
          res[[r]]$area    = area
          res[[r]]$ttarea  = ttarea
          res[[r]]$diff    = diff
          res[[r]]$maxdiff = maxdiff
          if(sortBy=="area")     res[[r]]=res[[r]][order(-area),]
          if(sortBy=="ttarea")   res[[r]]=res[[r]][order(-ttarea),]
          if(sortBy=="avg.diff") res[[r]]=res[[r]][order(-abs(diff)),]
          if(sortBy=="max.diff") res[[r]]=res[[r]][order(-abs(maxdiff)),]
          if(!is.null(removeIf)) res[[r]]=subset(res[[r]],subset=!eval(removeIf))
      } else {
           res[[r]]=data.frame(chr=character(0), start=numeric(0), end=numeric(0),
                               p1=numeric(0), p2=numeric(0), regionName=character(0),
                               indexStart=numeric(0), indexEnd=numeric(0), 
                               nprobes=numeric(0), area=numeric(0), ttarea=numeric(0), 
                               diff=numeric(0), maxdiff=numeric(0))
           if(is.null(p)) {
              colnames(res[[r]]) <- sub("p1", "m1", colnames(res[[r]]))
	      colnames(res[[r]]) <- sub("p2", "m2", colnames(res[[r]]))
           }
      }
      if(verbose) message(nrow(res[[r]])," DMR candidates found using cutoff=",cutoff[r],".")
  }
  if(verbose) message("\nDone")
  if(!paired) {
      return(list(tabs=res, p=p, l=l, chr=chr, pos=pos, pns=pns, 
              index=index, gm=lm, #controlIndex=controlIndex,
              groups=groups, args=args, comps=COMPS, package=package))
  } else {
      return(list(tabs=res, p=p, l=l, chr=chr, pos=pos, pns=pns, 
              index=index, DD=DD, sMD=sMD, #controlIndex=controlIndex,
              groups=groups, args=args, comps=COMPS, comps.names=COMPS.names, package=package))
  }
} # }}}

dmrFdr <- function(dmr, compare=1, numPerms=1000, seed=NULL, verbose=TRUE) {#{{{
	if (length(compare)!=1) stop("You must choose one comparison at a time when calculating FDRs. Please set dmr to be one of: ", 
	paste(names(dmr$tabs), collapse=", "), "\n")
	if (is.numeric(compare)) compare <- names(dmr$tabs)[compare]
	message("Calculating q-values for DMRs between ", compare, "\n")
	# Get probe order from TilingFeatureSet object
	pdInfo=get(dmr$package)
	class(pdInfo)="TilingFeatureSet" # Trick oligo so that pmChr, pmPosition work
	chr=pmChr(pdInfo)
	pos=pmPosition(pdInfo)
	chrpos <- paste(chr, pos)
	dchrpos <- paste(dmr$chr, dmr$pos)
	#idx <- which(!(chrpos %in% dchrpos))
	mis <- setdiff(chrpos, dchrpos)
	o <- order(chr, pos)
	chrpos <- chrpos[o]
	keepProbes <- !(chrpos %in% mis)
	o <- o[keepProbes]
	keep <- dmr$groups %in% unlist(strsplit(compare, "-"))
	# Recreate p or l with same sort order as in annotation package 
	if (is.null(dmr$p)) {
		l <- matrix(NA, nrow=length(pos), ncol=ncol(dmr$l))
		l[o,] <- dmr$l
		l <- l[,keep]
	} else {
		p <- matrix(NA, nrow=length(pos), ncol=ncol(dmr$p))
		p[o,] <- dmr$p
		p <- p[,keep]
	}
	n <- sum(keep)
	n1 <- sum(dmr$groups==unlist(strsplit(compare, "-"))[1])
	maxPerms <- choose(n, n1)
	if (numPerms=="all") numPerms <- maxPerms
	if (numPerms>maxPerms) {
		message("Given the sample sizes in the two groups the maximum number of permutations is ", maxPerms, sep="")
		numPerms <- maxPerms
	} 

	## Reshuffled group label DMRs
	if (!is.null(seed)) set.seed(seed)
	if (maxPerms<1e6) { # Enumerate all the combinations
		s <- sample(1:maxPerms, numPerms)
		grp1 <- combinations(n,n1)[s,]
	} else { 
		grp1 <- t(sapply(1:numPerms, function(x) sample(n,n1)))
	}

	if (verbose) message("Finding permuted data DMRs. Estimating time remaining")
	areas <- lapply(1:numPerms, function(i) {
		groups <- rep("grp2", n)
		groups[grp1[i,]] <- "grp1"
		if (is.null(dmr$p)) {
			st <- system.time(dmrPerm <- dmrFinder(dmr$package, l=l, 
				groups=groups, cutoff=dmr$args$cutoff, 
				filter=dmr$args$filter, ws=dmr$args$ws, 
				verbose=FALSE))[3]
		} else {
			st <- system.time(dmrPerm <- dmrFinder(dmr$package, p=p, 
				groups=groups, cutoff=dmr$args$cutoff, 
				filter=dmr$args$filter, ws=dmr$args$ws,
				verbose=FALSE))[3]
		}
		if (verbose & (i %in% round(seq(1, numPerms, length.out=10)))) {
			message(i, "/", numPerms, " (", prettyTime((numPerms-i)*st), 
				" remaining)", sep="")
		}
		dmrPerm$tabs[[1]][,dmr$args$sortBy]
	})

	nullDist <- unlist(areas)
	fn <- ecdf(nullDist)
	pval <- 1-fn(dmr$tabs[[compare]][,dmr$args$sortBy])
	pi0<-pi0.est(pval)$p0
	qval<-qvalue.cal(pval, pi0)
	dmr$tabs[[compare]] <- cbind(dmr$tabs[[compare]], pval, qval)
	if (!("numPerms" %in% names(dmr))) {
		dmr$numPerms <- rep(NA, length(dmr$tabs))
		names(dmr$numPerms) <- names(dmr$tabs)
	}
	dmr$numPerms[compare] <- numPerms
	return(dmr)
} # }}}

##The next 2 functions are for the paired=TRUE option of dmrFinder:

get.DD <- function(l, groups, compare, pairs){ # {{{

  ##Define COMPS the same as in get.tog():
  gIndex=split(seq(along=groups),groups)
  gIndex=gIndex[which(names(gIndex)%in%compare)]
  NAMES = names(gIndex)
  nums  <- match(compare,NAMES)
  COMPS <- matrix(nums,ncol=2,byrow=TRUE)
  COMPS.names <- t(apply(COMPS,1,function(x) NAMES[x]))
   
  DD    <- list() #unlike an array a list will accomodate drop-out
  for(i in 1:nrow(COMPS)){
    #Make DD[[i]], the differences for comparison i
    j=COMPS[i,1]
    k=COMPS[i,2]
    pd.j <- which(groups==NAMES[j])
    pd.k <- which(groups==NAMES[k])
    ord  <- match(pairs[pd.j],pairs[pd.k])
    if(all(is.na(ord))){
        message(paste("Comparison",i,"has no pairs! Ignoring this comparison."))
        DD[[i]] = NA
        next
    }
    if(any(is.na(ord))){
      pd.j <- pd.j[-which(is.na(ord))]
      ord  <- ord[-which(is.na(ord))]
    } #If any in j not in k, this will remove that one in j.
      #If any in k not in j, it's won't be included in pd.k anyway. 
    pd.k <- pd.k[ord]
    stopifnot(identical(pairs[pd.j],pairs[pd.k]))

    DD[[i]] <- as.matrix(l[,pd.j] - l[,pd.k])
    #if(absdiff) DD[[i]] <- abs(DD[[i]])^1
        #^.5 is down ladder of powers and makes symmetrical, but ^1 seems to work better.
        #log makes small differences large in magnitude (in the 
        #negative direction) which is not what we want.
    if(ncol(DD[[i]])==1) colnames(DD[[i]]) = colnames(l)[pd.j] #since no colname given
    names(DD)[i] = paste(NAMES[j],"-",NAMES[k],sep="")
  }
  return(list(DD=DD,COMPS=COMPS,COMPS.names=COMPS.names))
} # }}}

get.tt.paired <- function(DD,Indexes,filter=NULL,ws,verbose=TRUE) { # {{{
  #require(genefilter)
  ok = which(sapply(DD,length)>1)
  MD   <- matrix(NA,nrow=nrow(DD[[1]]),ncol=length(DD)) #mean differences
  vars <- matrix(NA,nrow=nrow(DD[[1]]),ncol=length(DD)) #their variances
  for(i in ok) {
    if(ncol(DD[[i]])<3) message(paste("Comparison",i,"is just using raw differences!  May want to set lower CUTOFF."))
    #if(ncol(DD[[i]])<3) stop(paste("Comparison",i,"has <3 pairs."))
    MD[,i]   <- rowMeans(DD[[i]])
    vars[,i] <- rowVars(DD[[i]])/ncol(DD[[i]]) #all NA if only 1 pair
  }

  #Now get smoothed MD and its vars:
  sMD <- matrix(NA,nrow=nrow(MD),ncol=ncol(MD))
  ses <- matrix(NA,nrow=nrow(vars),ncol=ncol(vars)) #myfilterse for new Charm returns sqrt(var)

  if(is.null(filter)) {
      Tukey = function(x) pmax(1 - x^2,0)^2
      fs= Tukey(seq(-ws,ws)/(ws+1));fs=fs/sum(fs)
  } else fs =filter

  if(verbose) message("Smoothing mean differences using window size of ",ws*2+1,":")
  if(verbose) pb = txtProgressBar(min=1, max=length(Indexes), initial=0, style=1)
  for(i in seq(along=Indexes)) {
    Index=Indexes[[i]]
    for(j in ok){
      sMD[Index,j] = myfilter( MD[Index,j],fs )
      if(ncol(DD[[j]])>2) ses[Index,j] = myfilterse( vars[Index,j],fs )
      if(ncol(DD[[j]])<3) ses[Index,j] = 1
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)
  colnames(MD)  = names(DD)
  colnames(sMD) = names(DD)
  #if(permute==FALSE & permute.paired==FALSE){
  #    save(MD,sMD,file=file.path(outpath,"MD.rda"),compress=TRUE) #for results.R
  #}
  #Finally make sT:
  sT <- sMD/ses
  #if(absdiff) sT = (sT - mean(sT))*(sT>mean(sT))
  return(list(sT=sT, MD=MD, sMD=sMD))
} # }}}

