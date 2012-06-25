###############################################################################
# Author: hui, huis@usc.edu
# Institute: USC, Epigenome Center
# Date creation: May 23, 2011
# Project Name: Hui_PhD_Analysis
###############################################################################
# Edits: tim, tim.triche@usc.edu
# What:  move lesions matrix to bottom
#        switch to WashU color scheme
#        other alterations?
###############################################################################

# one indicator
getBar <- function(SE, name, col1='gray', col2='black', col3=NULL) { # {{{
  if(is.logical(colData(SE)[[name]])) ifelse(colData(SE)[[name]], col2, col1)
  else ifelse(colData(SE)[[name]] %in% c('',' '), col1, col2)
} # }}}

# multiple indicators
getMatrix <- function(SE, names, col1='gray', col2='black', col3=NULL) { # {{{
  mat = do.call(rbind, lapply(names, function(x) getBar(SE, x, col1,col2,col3)))
  colnames(mat) <- colnames(SE)
  rownames(mat) <- names
  return(mat)
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
                        col = jet.colors(75),
                        ColSideColors, 
                        ColSideColors2,
                        RowSideColors,
                        rowColSrt=90,
                        colLabSrt=90,
                        cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc), 
                        labCex=1.5,
                        mtextCex=1.3,
                        ratio=1,
                        rowRatio=0.1,
                        colRatio=0.1,
                        labRow = NULL, 
                        labCol = NULL, 
                        main = NULL, 
                        xlab = NULL, 
                        ylab = NULL, 
                        keep.dendro = FALSE, 
                        Rdend=TRUE,
                        Cdend=TRUE,
		                    verbose = getOption("verbose"),
                        output=FALSE, 
                        ...) {

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

	if (is.null(Rowv)) 
		Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) 
		Colv <- colMeans(x, na.rm = na.rm)
		
	if (!is.null(order.x)){
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
		} else{
		if (doCdend) {
			if (inherits(Colv, "dendrogram")) 
				ddc <- Colv
			else if (identical(Colv, "Rowv")) {
				if (nr != nc) 
					stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
				ddc <- ddr
			}else {
				hcc <- hclustfun(distfun(if (symm) 
											x
										else t(x)))
				ddc <- as.dendrogram(hcc)
				if (!is.logical(Colv) || Colv) 
					ddc <- reorderfun(ddc, Colv)
			}
			if (nc != length(colInd <- order.dendrogram(ddc))) 
				stop("column dendrogram ordering gave index of wrong length")
		}
		else colInd <- 1:nc
	}
	
	if (!is.null(order.y)){
		doRdend<-FALSE
		IDs<-split(1:nr,order.y)
		rowInd<-NULL
		for (i in 1:length(table(order.y))){
			ID<-IDs[[i]]
			y.i<-x[ID,]
			hcr.i <- hclustfun(distfun(y.i))
			ddr.i <- as.dendrogram(hcr.i)
			rowInd<-c(rowInd,ID[order.dendrogram(ddr.i)])
		}
	}else{
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
		
	} 
	
	

	
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) 
				if (is.null(rownames(x))) 
					(1:nr)[rowInd]
				else rownames(x)
			else labRow[rowInd]
	labCol <- if (is.null(labCol)) 
				if (is.null(colnames(x))) 
					(1:nc)[colInd]
				else colnames(x)
			else labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
		sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
		sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	lmat <- rbind(c(NA, 3), c(2:1))
	lwid <- c(if (Rdend) 1 else 0.25, 4)
	lhei <- c((if (Cdend) 1 else 0.05), 
			4)
	if (!missing(ColSideColors)) {
		if (!is.matrix(ColSideColors)) 
			stop("'ColSideColors' must be a matrix")
		if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
				nc) 
			stop("'ColSideColors' dim()[2] must be of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1], 4*rowRatio, lhei[2])
	}
	if (!missing(RowSideColors)) {
		if (!is.matrix(RowSideColors)) 
			stop("'RowSideColors' must be a matrix")
		if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
				nr) 
			stop("'RowSideColors' must be a character vector of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1], 4*colRatio, lwid[2])
	}
	if(!is.null(x2)){
		lmat <- cbind(lmat,c(rep(NA,nrow(lmat)-1),max(lmat,na.rm=T)+1))
		lwid <- c(lwid,4*ratio)
		
	}
	lmat[is.na(lmat)] <- 0
		if (!missing(ColSideColors2)) {
		lmat [which.max(lmat)-1]<- max(lmat)+1
			}
	
	
	
	if (verbose) {
		cat("layout: widths = ", lwid, ", heights = ", lhei, 
				"; lmat=\n")
		print(lmat)
	}
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	if (!missing(RowSideColors)) {
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
		image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
		if (length(colnames(RowSideColors)) > 0) {
			axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), labels=FALSE, las = 2, tick = FALSE)
		text(0:(dim(rsc)[2]-1)/(dim(rsc)[2] - 1) , -0.02, colnames(RowSideColors), 
		
					srt = rowColSrt,  adj=c(1,0.5), xpd=TRUE, font=2,cex=labCex)
			if (!is.null(mtextRow)){
				mtext(mtextRow,side=2,line=1,font=2,cex=mtextCex)
				
			}
			
		}
	}
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
		image(csc, col = as.vector(csc.colors), axes = FALSE)
		if (length(colnames(ColSideColors)) > 0) {
			
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
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), col=col, axes = FALSE, xlab = "", ylab = "")
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
	if (!is.null(x2)){
		par(mar=c(margins[1],5,sidemar,margins[2]))
		rownames.ordered<-row.x[rowInd]
		newrow<-rownames.ordered
		hcc2 <- hclustfun(distfun(if (symm) 
									x2
								else t(x2)))
		ddc2 <- as.dendrogram(hcc2)
		newcol<-order.dendrogram(ddc2)
		
		x2.ordered<-as.matrix(x2[newrow,newcol])
		nc2<-dim(x2)[2]
		image(1:nc2, 1:nr, t(x2.ordered), xlim =0.5+c(0, nc2), ylim = 0.5+c(0, nr), axes = FALSE,xlab = "", ylab = "",col=col,...)
		#axis(1, 1:nc2, labels = FALSE, las = 2, line = -0.5, tick = 0) #, cex.axis = cexCol)
	#text(c(1:nc2), 1, labels = colnames(x2)[newcol], srt = colLabSrt, adj=c(1.5,1.5), xpd=TRUE, offset=2, font=2,cex=cexCol)
		if (!is.null(mtext2)){
		mtext(mtext2,1, adj=0.5, line=mline2,font=2,cex=mtextCex)
	}
	}
	
	if (!missing(ColSideColors2)) {
		par(mar = c(0.5,5, 0, margins[2]))
		csc = ColSideColors2[newcol, ]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names) {
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(csc, col = as.vector(csc.colors), axes = FALSE)
		
	}
	
	
	
	if (output){
		
	out<-list()
	out$colInd<-colInd
	out$rowInd<-rowInd
	out$colInd2<-newcol
	out$rowInd2<-newrow
return(out)
	}
}



############################################################################
# HISTORY:
# May 23, 2011
# o Created. by hui
############################################################################
