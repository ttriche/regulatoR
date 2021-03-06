###############################################################################
# Author: hui, huis@usc.edu
# Institute: USC, Epigenome Center
# Date creation: May 23, 2011
# Project Name: Hui_PhD_Analysis
###############################################################################
# Edits: tim, tim.triche@usc.edu
###############################################################################

# standardized mutation matrix colors
std.cols <- c( # {{{
              '#E0E0E0', # gray
              '#000000', # black
              '#BB0000', # red
              '#777777', # dark gray
              '#008800', # green
              '#888800', # yellow
              '#008888', # teal
              '#880088', # purple
              '#0000AA', # blue
              '#004488', # tealish
              '#886644', # tealish
              '#FFFF00', # tealish
              '#0088FF', # tealish
              '#664422', # tealish
              '#FF00FF' # tealish
              ) # }}}
alt.cols <- c(  # {{{
              '#E0E0E0', # gray
              '#000088', # blue  (fused)
              '#000000', # black (mutated)
              '#008800', # green
              '#888800', # yellow
              '#008888', # teal
              '#880088'  # purple
              ) # }}}

# standardized fonts 
type.font <- c('cell' = 1,
               'fusion' = 1,
               'gene' = 3)

std.anno <- c('CD34.1',
              'CD34.2',
              'CD34.3')
new.anno <- c(std.anno,
              'PROS.4',
              'PROS.5',
              'PROS.6',
              'PMN.4',
              'PMN.5',
              'PMN.6',
              'MONO.4',
              'MONO.5',
              'MONO.6')
#std.anno <- c(std.anno, 'SPACER', 'CPGI')
std.anno <- c(new.anno, 'SPACER', 'CPGI')
rev.anno <- rev(std.anno) # for low-CpG panel

# normal bone marrow covariates
ctls.type <- c('PMN' = 'cell',
               'CD34' = 'cell',
               'CD14' = 'cell',
               'CD15' = 'cell')
ctls.std <- names(ctls.type)

# LAML mutations/fusions of note
# formatting is ala Tim Ley
old.muts.type <- c( # {{{
               #'kEquals3' = 'cluster',
               #'kEquals4' = 'cluster',
               #'kEquals5' = 'cluster',
               #'kEquals6' = 'cluster',
               #'kEquals11' = 'cluster',
               #'cluster' = 'cluster',
               # 'FAB' = 'cluster',

               'SPACER' = 'spacer',

               'DNMT3A' = 'gene',
               'NPM1' = 'gene',
               'FLT3' = 'gene',
               #'NRAS.KRAS' = 'gene',
               #'Phosphatase' = 'gene', ## PTPN11, PTPRT, PTPRN, PTPRD
                                       ### PTPRJ, PTPRG, PTPRE, PTPRN2
               #'ABL.CBL.KIT' = 'gene',
               'TET1.TET2' = 'gene',
               'IDH1.IDH2' = 'gene', 
               #'Polycomb' = 'gene', ## EZH2, EED, SUZ12, ASXL1, CBX7
               #'ETV6.IKZF1.IKZF4' = 'gene',
               #'Spliceosome' = 'gene', ## SF3B1, PRPF8, PRPF4B, U2AF1, U2AF2,
                                       ### SRSF6, TRA2B, CSTF2T, DDX1, DDX23, 
                                       ### DHX32, METTL3, PLRG1, PRPF3, RBMX,
                                       ### SNRNP200, SRRM1, SRRM2, SUPT5H,
                                       ### U2AF1L4
               #'Cohesin' = 'gene', ## STAG2, RAD21, SMC3, SMC1A
               'RUNX1' = 'gene', 
               'CEBPA' = 'gene', 
               'TP53' = 'gene',

               'SPACER' = 'spacer',
                
               'MLL.fusions' = 'fusion', ## several

               'SPACER' = 'spacer',

               'RUNX1.RUNX1T1' = 'fusion',
               'MYH11.CBFB' = 'fusion',
               'PML.RARA' = 'fusion',
               'PICALM.MLLT10' = 'fusion' ## two
               ) # }}}
muts.type <- rev( # {{{
                 c(
                   M0 = 'gene',
                   M1 = 'gene',
                   M2 = 'gene',
                   M3 = 'gene',
                   M4 = 'gene',
                   M5 = 'gene',
                   M6 = 'gene',
                   M7 = 'gene',
                   SPACER = 'spacer',
                   PML.RARA = 'fusion',
                   MYH11.CBFB = 'fusion',
                   RUNX1.RUNX1T1 = 'fusion',
                   PICALM.MLLT10 = 'fusion',
                   SPACER = 'spacer',
                   NPM1 = 'gene',
                   TP53 = 'gene',
                   WT1 = 'gene',
                   DNMT3A = 'gene',
                   DNMT3B = 'gene',
                   DNMT1 = 'gene',
                   TET1 = 'gene',
                   TET2 = 'gene',
                   IDH1 = 'gene',
                   IDH2 = 'gene',
                   FLT3 = 'gene',
                   KIT = 'gene',
                   Other_Tyr_kinases = 'gene',
                   Ser.Thr_kinases = 'gene',
                   KRAS.NRAS = 'gene',
                   PTPs = 'gene',
                   RUNX1 = 'gene',
                   CEBPA = 'gene',
                   Other_myeloid_TFs = 'gene',
                   MLL_fusions = 'fusion',
                   MLL_PTD = 'gene',
                   NUP98.NSD1 = 'fusion',
                   ASXL1 = 'gene',
                   EZH2 = 'gene',
                   KDM6A = 'gene',
                   Other_modifiers = 'gene',
                   Cohesin = 'gene',
                   Spliceosome = 'gene',
                   SPACER = 'spacer',
                   Cytogenetic_risk = 'risk'
                )
              ) # }}}
# muts.std <- names(rev(muts.type))
muts.std <- c('SPACER'='SPACER','SPACER'='SPACER')
muts.type <- muts.std ## doh

# traditional
jet.colors <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F", 
                                 "yellow", "#FF7F00", "red", "#7F0000"))

# monochrome
ley.colors <- colorRampPalette(c("white","black"))

# Expression
exp.colors <- colorRampPalette(c("green","black","red"))

# Age
age.colors <- colorRampPalette(c("green","yellow","red"))

# Risk
risk.colors = function(...) return(c("white","pink","red","gray50"))

# accessibility
dgf.colors <- colorRampPalette(
                c(rep("white",6),"yellow","orange","red","darkred")) 

# accessibility
dhs.colors <- colorRampPalette(c('white','darkred'))

# one indicator, getBar(SE, 'covariate') OR getBar(SE$covariate) will work too
getBar <- function(SE, name=NULL, cols=std.cols, alt=alt.cols) { # {{{
  if( toupper(name) == 'SPACER' ) return(c(rep('white', ncol(SE))))
  if( options()$verbose == TRUE ) message(name)
  if(is(SE, 'SummarizedExperiment')) {
    x <- colData(SE)[[name]]
  } else {
    x <- SE ## usually will be SE$something instead
  }
  if(toupper(name)=='GENDER') {
    ifelse(colData(SE)[[name]]=='M','lightblue','pink')
  } else if(toupper(name)=='CYTOGENETIC_RISK') {
    risk.colors()[colData(SE)[[name]]]
  } else if(toupper(name)=='AGE') {
    age.colors(100)[colData(SE)[[name]]]
  } else if(is.logical(x)) {
    ifelse(colData(SE)[[name]], cols[2], cols[1])
  } else if(is.factor(x)) {
    cols[as.numeric(colData(SE)[[name]])]
  } else {
    ifelse(x %in% c('',' ','false','FALSE'), cols[1], cols[2])
  }
} # }}}
getMatrix <- function(SE, muts, cols=std.cols, alt=alt.cols) { # {{{
  mat <- matrix(NA, ncol=ncol(SE), nrow=length(muts))
  colnames(mat) <- colnames(SE)
  rownames(mat) <- names(muts)
  for(i in seq_along(muts)) {
    coloring = cols
    x = names(muts)[i]
    message(x)
    if(muts[i] == 'fusion') {
      coloring = alt
      mat[ x, ] = getBar(SE, x, coloring)
      rownames(mat)[i] <- gsub('\\.', '-', x)
    } else if(muts[i] == 'gene') {
      mat[ x, ] = getBar(SE, x, coloring)
      rownames(mat)[i] <- gsub('\\.','/', x)
    } else if(muts[i] == 'spacer') {
      mat[ x, ] = getBar(SE, x, c('white','white'))
    } else {
      mat[ x, ] = getBar(SE, x, coloring)
    }
    rownames(mat)[i] <- gsub('_',' ', rownames(mat)[i])
    rownames(mat)[i] <- gsub('SPACER',' ', rownames(mat)[i])
  }
  return(mat)
} # }}}

getProbeBar <- function(SE, name=NULL, col1='#9F9FA3',col2='#000000') { # {{{
  if(is(SE, 'SummarizedExperiment')) x <- values(SE)[[name]]
  else x <- SE ## usually will be values(rowData(SE))$something instead
  if( substr(toupper(name), 1, 6) == 'SPACER' ) {
    return( rep('white', nrow(SE)) )
  } else if(substr(toupper(name), 1, 3)=='CPG') {
    ifelse(x=='island','darkgreen',ifelse(x=='shore','green','white'))
  } else if(toupper(name) %in% c('OECG','DGF','CD34.DGF','LOG.DGF','LADS')) {
    dgf.colors(100)[x]
  } else if(toupper(name) %in% c('CD34.DHS','CD14.DHS','K562.DHS',
                                 'CMK.DHS','NB4.DHS','HL60.DHS','DHS')) {
    ley.colors(100)[x / 2]
  } else if(toupper(name) %in% c('CD34.H3K4','CD34.H3K27','CD34.H3K36',
                                 'H3K4ME3','H3K27ME3','H3K36ME3',
                                 'H3K4ME1','H3K27AC')) {
    ley.colors(100)[x * 10]
  } else if(toupper(name) %in% c('CD34.1','CD34.2','CD34.3','CD133','PBMC',
                                 'K562','CMK','NB4','HL60','MONO.4','MONO.5',
                                 'MONO.6','PROS.4','PROS.5','PROS.6','PMN.1',
                                 'PMN.2','PMN.3','PMN.4','PMN.5','PMN.6')) {
    jet.colors(100)[x * 100]
  } else {
    ifelse(values(rowData(SE))[[name]], col2, col1)
  }
} # }}}
getProbeMatrix <- function(SE, covs) { # {{{
  mat <- matrix( rep('#FFFFFF', nrow(SE)), ncol=1 )
  mat = cbind(mat,
              do.call(cbind, lapply(covs, function(x) getProbeBar(SE, x))))
  colnames(mat) <- c('', gsub('CD34.DGF','DGF', gsub('SPACER.*',' ',covs)))
  colnames(mat) <- gsub('CD34.*','CD34', colnames(mat))
  colnames(mat) <- gsub('MONO.*','MONO', colnames(mat))
  colnames(mat) <- gsub('PROS.*','PROS', colnames(mat))
  colnames(mat) <- gsub('PMN.*','PMN', colnames(mat))
  rownames(mat) <- rownames(SE)
  return(mat)
} # }}}

rowSdSdMax <- function(x, assay=NULL) { # {{{
  require(matrixStats)
  if(is(x, 'SummarizedExperiment')) x <- asy.fast(x, assay)
  return(rowSds(x,na.rm=T)/sqrt(rowMeans(x,na.rm=T)/(1-(rowMeans(x, na.rm=T)))))
} # }}}

by.sd <- function(x, howmany=1000, assay=NULL) { # {{{
  if('sd' %in% names(values(x))) {
    return(x[ head(order(values(rowData(x))$sd, decreasing=TRUE), howmany), ])
  } else {
    require(matrixStats)
    values(x)$sd <- rowSds(asy.fast(x, assay), na.rm=TRUE)
    return(x[ head(order(values(rowData(x))$sd, decreasing=TRUE), howmany), ])
  }
} # }}}

by.sdSdMax <- function(SE, howmany=1000, assay=NULL) { # {{{
  if(!('sdSdMax' %in% names(values(rowData(SE))))) {
    values(rowData(SE))$sdSdMax <- rowSdSdMax(asy.fast(SE, assay))
  }
  return(SE[head(order(values(rowData(SE))$sdSdMax,decreasing=TRUE), howmany),])
} # }}}

by.mad <- function(SE, howmany=1000, assay=NULL) { # {{{
  if(!('mad' %in% names(values(rowData(SE))))) {
    require(matrixStats)
    values(rowData(SE))$mad <- rowMads(asy.fast(SE, assay), na.rm=TRUE)
  }
  return(SE[ head(order(values(rowData(SE))$mad, decreasing=TRUE), howmany), ])
} # }}}

fetchby <- function(SE, how, howmany, CpH=FALSE) { # {{{
  if(any(grepl('^ch',rownames(SE)))) {
    if(CpH==TRUE) se <- split(SE, substr(rownames(SE),1,2))$ch
    else se <- split(SE, substr(rownames(SE),1,2))$cg
  } else if(CpH==TRUE) {
    stop('No CpH probes were found.')
  } else if(any(grepl('^rs', rownames(SE)))) {
    se <- split(SE, substr(rownames(SE),1,2))$cg
  } else {
    se <- SE
  }
  xx <- switch(how, sd=by.sd(se, howmany), 
                    mad=by.mad(se, howmany), 
                    sdmax=by.sdSdMax(se, howmany))
  x <- asy.fast(xx)
  rownames(x) <- rownames(xx)
  colnames(x) <- colnames(xx)
  rm(xx)
  gc()
  return(x)
} # }}}

# how lazy am I? VERY.
coolmap <- function(SE1,
                    SE2=NULL,
                    how='sd',
                    howmany=1000,
                    method='ward',
                    muts=muts.type, 
                    ctls=ctls.type,
                    col=jet.colors(255),
                    logit=FALSE,
                    Rdend=FALSE,
                    Cdend=FALSE,
                    probeAnno=std.anno,
                    output=TRUE,
                    labRow='',
                    CpH=FALSE,
                    k=15,
                    rowRatio=0.2,
                    colRatio=0.2,
                    ...) 
{ # {{{

  exclude <- c()
  if(any(seqnames(SE1) %in% c('chrX','chrY'))) { # {{{ tidy up
    exclude <- union(exclude, which(seqnames(SE1) %in% c('chrX','chrY')))
  } # }}}
  if(any(rowSums(is.na(asy.fast(SE1)))/ncol(SE1) >= 0.5)) { # {{{
    exclude = union(exclude,which(rowSums(is.na(asy.fast(SE1)))/ncol(SE1)>=.5))
  } # }}}

  ## this needs revamping due to the repeat situation
  if(!exists('HM450.SNP15.HG19')) data(HM450.SNP15.HG19)
  exclude <- union(exclude, which(rowData(SE1) %in% HM450.SNP15.HG19))
  if(!exists('HM450.RPT15.HG19')) data(HM450.RPT15.HG19)
  exclude <- union(exclude, which(rownames(SE1) %in% names(HM450.RPT15.HG19)))
  if('sd' %in% names(values(SE1))) {
    exclude <- union(exclude, which(is.na(values(SE1)$sd)))
  }
  require(impute)
  x <- x2 <- NULL
  x <- fetchby(SE1[ -exclude, ], how, howmany, CpH=CpH)
  if(anyMissing(x)) x <- impute.knn(x)$data

  if(!is.null(probeAnno)) {
    RowSideColors = getProbeMatrix(SE1[rownames(x),], probeAnno)
  } else {
    RowSideColors = NULL
  }

  if(!is.null(SE2)) {
    x2 <- asy.fast(SE2[ rownames(x), ])
    rownames(x2) <- rownames(x)
    colnames(x2) <- colnames(SE2)
    if(anyMissing(x2)) x2 <- impute.knn(x2)$data
  }

  z <- z2 <- NULL
  if(!is.null(muts)) z <- t(getMatrix(SE1, muts))
  if(!is.null(SE2) && !is.null(ctls)) z2 <- t(getMatrix(SE2, ctls))

  hf <- function(w) hclust(w, method=method)
  
  if(logit) x <- logit(x)
  if(logit && !is.null(x2)) x2 <- logit(x2)
  
  if(is.null(SE2)) {
    labCol = '' #  gsub('TCGA(_|-AB-)','',colnames(SE1))
    if(is.null(z)) {
      out <- heatmap.minus(x=x, col=col, hclustfun=hf, scale='none', 
                           Rdend=Rdend, Cdend=Cdend, labRow=labRow, 
                           labCol=labCol, RowSideColors=RowSideColors, 
                           x2names=FALSE, output=TRUE, rowRatio=rowRatio,
                           colRatio=colRatio, k=k, ...)
    } else {
      out <- heatmap.minus(x=x, ColSideColors=z, col=col, hclustfun=hf, 
                           scale='none', Rdend=Rdend, Cdend=Cdend, 
                           labRow=labRow, labCol=labCol, 
                           RowSideColors=RowSideColors, 
                           x2names=FALSE, output=TRUE, rowRatio=rowRatio,
                           colRatio=colRatio, k=k, ...)

    }
  } else { 
    out <- heatmap.minus(x=x, x2=x2, ColSideColors=z, ColSideColors2=z2, 
                         col=col, hclustfun=hf, scale='none', Rdend=Rdend, 
                         Cdend=Cdend, RowSideColors=RowSideColors, 
                         x2names=FALSE, labRow=labRow, output=TRUE, 
                         rowRatio=rowRatio, colRatio=colRatio, k=k, ...)
  }
  
  out$colData <- colData(SE1)[out$colNames,]
  out$rowData <- rowData(SE1)[out$rowNames]
  if(Cdend==TRUE) rownames(out$clusters) <- colnames(SE1)
  return(out)

} # }}}

# x is tumor, x2 is normal
heatmap.minus<-function(x,
                        x2=NULL, 
                        order.x=NULL,
                        order.y=NULL,
                        Rowv=NULL, 
                        Colv = if (symm) "Rowv" else NULL, 
                        mtext1=NULL,
                        mtext2=NULL, 
                        mline1=2,
                        mline2=2, 
                        mtextRow=NULL, 
                        distfun = dist, 
                        hclustfun = hclust, 
                        reorderfun = function(d, w) { reorder(d, w) }, 
                        add.expr, 
                        symm = FALSE, 
                        revC = identical(Colv, "Rowv"), 
                        scale = c("none", "row", "column"), 
                        na.rm = TRUE, 
                        margins = c(6, 6),
                        sidemar = 0.5, 
                        splitmar = 2,
                        col = jet.colors(75),
                        ColSideColors, 
                        ColSideColors2,
                        RowSideColors,
                        rowColSrt=90,
                        colLabSrt=90,
                        cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc), 
                        labCex=1.25,
                        mtextCex=1.25,
                        ratio=1,
                        rowRatio=0.1,
                        colRatio=0.1,
                        labRow = NULL, 
                        labCol = NULL, 
                        main = NULL, 
                        xlab = NULL, 
                        ylab = NULL, 
                        keep.dendro = FALSE, 
                        Rdend=FALSE,
                        Cdend=TRUE,
		                    verbose = getOption("verbose"),
                        output=FALSE, 
                        x2names=FALSE,
                        revT=FALSE,
                        k=15,
                        leftLabels=FALSE,
                        useRaster=FALSE,
                        ...) 
{ # {{{

	scale <- ifelse(symm && missing(scale), "none", match.arg(scale))
	if(length(di<-dim(x))!=2 || !is.numeric(x)) stop("x must be a numeric matrix")
  doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	nr <- di[1]
	nc <- di[2]
	row.x<-rownames(x)
	if (nr <= 1 || nc <= 1) 
		stop("'x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2) 
		stop("'margins' must be a numeric vector of length 2")

	if (is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)
		
	if (!is.null(order.x)){ # {{{
		Cdend<-FALSE
		doCdend<-FALSE
		IDs<-split(1:nc,order.x)
		colInd<-NULL
		for (i in 1:length(table(order.x))){
			ID<-IDs[[i]]
			x.i<-x[,ID]
			hcc.i <- hclustfun(distfun(t(x.i)))
			ddc.i <- as.dendrogram(hcc.i)
			colInd<-c(colInd,ID[order.dendrogram(ddc.i)])
		}
	} else {
		if (doCdend) {
			if (inherits(Colv, "dendrogram")) 
				ddc <- Colv
			else if (identical(Colv, "Rowv")) {
				if (nr != nc) 
					stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
				ddc <- ddr
			} else {

        ## automatic cluster calls here
        hcc <- hclustfun(distfun(if(symm) x	else t(x)))
        sizes <- seq(2, k, 1)
        names(sizes) <- paste0('kEquals', sizes)
        clusters <- do.call(cbind, lapply(sizes, function(i) cutree(hcc, k=i)))
        clusters <- DataFrame(clusters)
        for(kk in names(clusters)) clusters[, kk] <- as.factor(clusters[, kk])
        ## done with automatic cluster calls 

        ddc <- as.dendrogram(hcc)
				if(!is.logical(Colv) || Colv) 
					ddc <- reorderfun(ddc, Colv)
			}
			if (nc != length(colInd <- order.dendrogram(ddc))) 
				stop("column dendrogram ordering gave index of wrong length")
		}
		else colInd <- 1:nc
	} # }}}
	if (!is.null(order.y)){ # {{{
		doRdend<-FALSE
		IDs<-split(1:nr,order.y)
		rowInd<-NULL
		for (i in 1:length(table(order.y))){
			ID<-IDs[[i]]
			y.i<-x[ID,]
			hcr.i <- hclustfun(distfun(y.i))
			ddr.i <- as.dendrogram(hcr.i)
			rowInd<-c(rowInd,ID[order.dendrogram(ddr.i)])
		} # }}}
	}else{ # {{{
		if (doRdend) {
			if (inherits(Rowv, "dendrogram")) 
				ddr <- Rowv
			else {
				hcr <- hclustfun(distfun(x))
				ddr <- as.dendrogram(hcr)
				if (!is.logical(Rowv) || Rowv) 
					ddr <- reorderfun(ddr, Rowv)
			}
			if (nr != length(rowInd <- order.dendrogram(ddr))) 
				stop("row dendrogram ordering gave index of wrong length")
		}
		else rowInd <- 1:nr
		
	}  # }}}
	
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) { # {{{
      				if (is.null(rownames(x))) {
                (1:nr)[rowInd]
      				} else {
                rownames(x)
              }
			      } else {
              labRow[rowInd]
            } # }}}
	labCol <- if (is.null(labCol)) { # {{{
      				if (is.null(colnames(x))) {
      					(1:nc)[colInd]
      				} else {
                colnames(x)
			        }
            } else {
              labCol[colInd]
            } # }}}
	if (scale == "row") { # {{{
		x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
		sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	} # }}}
	else if (scale == "column") { # {{{
		x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
		sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	} # }}}
	lmat <- rbind(c(NA, 3), c(2:1))
	lwid <- c(ifelse(Rdend, 1, 0.5), 4)
	lhei <- c(ifelse(Cdend, 1, 0.25), 4)
	if (!missing(ColSideColors)) { # {{{
		if (!is.matrix(ColSideColors)) 
			stop("'ColSideColors' must be a matrix")
		if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
				nc) 
			stop("'ColSideColors' dim()[2] must be of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1], 4*rowRatio, lhei[2])
	} # }}}
	if (!missing(RowSideColors)) { # {{{
		if (!is.matrix(RowSideColors)) 
			stop("'RowSideColors' must be a matrix")
		if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
				nr) 
			stop("'RowSideColors' must be a character vector of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1], 4*colRatio, lwid[2])
	} # }}}
	if(!is.null(x2)){ # {{{
		lmat <- cbind(lmat,c(rep(NA,nrow(lmat)-1),max(lmat,na.rm=T)+1))
		lwid <- c(lwid,4*ratio)
	} # }}}
	lmat[is.na(lmat)] <- 0
  if (!missing(ColSideColors2)) lmat[which.max(lmat)-1] <- max(lmat)+1
	if (verbose) { # {{{
		cat("layout: widths = ", lwid, ", heights = ", lhei, 
				"; lmat=\n")
		print(lmat)
	} # }}}
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	if (!missing(RowSideColors)) { # {{{
		par(mar = c(margins[1], 2, sidemar, 0))
		rsc = RowSideColors[rowInd, ]
		rsc.colors = matrix()
		rsc.names = names(table(rsc))
		rsc.i = 1
		for (rsc.name in rsc.names) {
			rsc.colors[rsc.i] = rsc.name
			rsc[rsc == rsc.name] = rsc.i
			rsc.i = rsc.i + 1
		}
		rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
		image(t(rsc), col = as.vector(rsc.colors), axes=FALSE, useRaster=useRaster)
		if (length(colnames(RowSideColors)) > 0 && leftLabels == TRUE ) {
			axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), labels=FALSE, las = 2, tick = FALSE)
		text(0:(dim(rsc)[2]-1)/(dim(rsc)[2] - 1) , -0.02, colnames(RowSideColors), 
		
					srt = rowColSrt,  adj=c(1,0.5), xpd=TRUE, font=2,cex=labCex)
			if (!is.null(mtextRow)){
				mtext(mtextRow,side=2,line=1,font=2,cex=mtextCex)
			}
		}
	} # }}}
	if (!missing(ColSideColors)) {
		par(mar = c(0.5, sidemar, 0, if (is.null(x2)) margins[2] else 0))
		csc = ColSideColors[colInd, ]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names) {
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(csc, col = as.vector(csc.colors), axes = FALSE, useRaster=useRaster)
		if (length(colnames(ColSideColors)) > 0 && leftLabels == TRUE) {
			
			axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors), font=2, las = 2, tick = FALSE,cex.axis=labCex)
		}
	}
	par(mar = c(margins[1], sidemar, sidemar, if (is.null(x2)) margins[2] else 0))
	if (!symm || scale != "none") {
		x <- t(x)
	}
	if (revC) {
		iy <- nr:1
		ddr <- rev(ddr)
		x <- x[, iy]
	}
	else iy <- 1:nr
  if(revT==TRUE) {
    image(1:nc, 1:nr, x[nc:1,], xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), 
          col=col, axes = FALSE, xlab = "", ylab = "", useRaster=useRaster)
  } else { 
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), 
          col=col, axes = FALSE, xlab = "", ylab = "", useRaster=useRaster)
  } 
  axis(1, 1:nc, labels = FALSE, las = 2, line = -0.5, tick = 0) #, cex.axis = cexCol)
	text(c(1:nc), 1, labels = labCol, srt = colLabSrt, adj=c(1.5,1.5), xpd=TRUE, offset=2, font=2,cex=cexCol)
	if (!is.null(mtext1)){
	mtext(mtext1,1,line=mline1,cex=mtextCex,font=2)
}
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexRow)
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	par(mar = c(margins[1], 0, sidemar, 0))
	if (Rdend) 
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	else frame()
	par(mar = c(0, sidemar, if (!is.null(main)) 5 else 0, if (is.null(x2)) margins[2] else 0))
	if (Cdend) 
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	else frame() 
		if (!is.null(main)) 
		frame()
	if (!is.null(main)) 
		title(main, cex.main = 1.5 * op[["cex.main"]])
	invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
							doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))

	if (!is.null(x2)){ # {{{
		par(mar=c(margins[1], splitmar, sidemar,margins[2]))
		rownames.ordered<-row.x[rowInd]
		newrow<-rownames.ordered
		hcc2 <- hclustfun(distfun(if (symm) 
									x2
								else t(x2)))
		ddc2 <- as.dendrogram(hcc2)
		newcol<-order.dendrogram(ddc2)
		
		x2.ordered<-as.matrix(x2[newrow,newcol])
		nc2<-dim(x2)[2]
		image(1:nc2, 1:nr, t(x2.ordered), xlim =0.5+c(0, nc2), ylim = 0.5+c(0, nr),
          axes = FALSE, xlab = "", ylab = "", col=col, useRaster=useRaster, ...)
	  axis(1, 1:nc2, labels = FALSE, las = 2, line = -0.5, tick = 0) #, cex.axis = cexCol)
	  if(x2names != FALSE) text(c(1:nc2), 1, labels = colnames(x2)[newcol], srt = colLabSrt, adj=c(1.5,1.5), xpd=TRUE, offset=2, font=2,cex=cexCol)
		if (!is.null(mtext2)){
	  	mtext(mtext2,1, adj=0.5, line=mline2,font=2,cex=mtextCex)
	  }
	} # }}}
	if (!missing(ColSideColors2) && !is.null(ColSideColors2)) { # {{{
		par(mar = c(0.5, splitmar, 0, margins[2]))
		csc = ColSideColors2[newcol, ]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names) {
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
    colnames(csc) <- colnames(ColSideColors2)
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(csc, col = as.vector(csc.colors), axes = FALSE)
    axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors2), font=2, las = 2, tick = FALSE,cex.axis=labCex)
	} # }}}
	
	if (output){
    out<-list()
    if(revT==TRUE) {
      out$colInd <- rev(colInd)
      out$colNames <- rev(rownames(x))
      out$clusters <- rev(clusters)
    } else { 
      out$colInd <- colInd
      out$colNames <- rownames(x)
      out$clusters <- clusters
    }
    out$ddc <- ddc
    out$rowInd <- rowInd
    out$rowNames <- colnames(x)
    out$ddr <- ddr
	  if (!missing(ColSideColors2) && !is.null(ColSideColors2)) {
      out$colInd2<-newcol
      out$rowInd2<-newrow
    }
    return(out)
	}
} # }}}
