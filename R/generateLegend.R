###############################################################################
# Author: hui, huis@usc.edu
# Institute: USC, Epigenome Center
# Date creation: Jun 7, 2011
# Project Name: Hui_PhD_Analysis
# TODO: TODO
###############################################################################

getPal<-function(palName="Paired",start=1,step=1,nullCol="white",n=8,na=FALSE){
	library("RColorBrewer")
	thisPal<-brewer.pal(n, palName)[c(start:n,(1:(start-1)))]
	
	if (na){
		
		thisPal<-c("white",thisPal)
	}
	
	if (step !=1){
		i<-1:2
		thisPal<-thisPal[1+(step)*(i-1)]
	}
	thisPal
}



getLegend<-function(labels,sublabels,continuous,subcols,cex.title=1, cex.sub=0.8) {
n<-length(labels)
n1<-unlist(lapply(sublabels,length))
nsubcol<-ceiling(n1/2)
n.add<-cumsum(c(0,nsubcol[-length(nsubcol)]))

y.unit<-1/(n+sum(nsubcol))
bord<-y.unit*0.618

plot(0:1,0:1,type="n",axes=FALSE,xlab="",ylab="")

for (i in 1:n){
	text(0,1-y.unit*(i-1+n.add[i]),toupper(labels[i]),font=2,xpd=TRUE,cex=cex.title,adj=c(0,0))
	if (continuous[i]=="nc"){
		
		for (j in seq(1,n1[i],2)){
			rect(0,1-y.unit*(i+n.add[i]+floor(j/2)),bord,1-y.unit*(i+n.add[i]+floor(j/2))+bord,col=subcols[[i]][j],border=TRUE)
			text(bord+0.01,1-y.unit*(i+n.add[i]+floor(j/2))+(y.unit-bord)/2,sublabels[[i]][j],cex=cex.sub,adj=c(0,0.5),font=2)
			if (j<n1[i]){
				rect(0.3,1-y.unit*(i+n.add[i]+floor(j/2)),0.3+bord,1-y.unit*(i+n.add[i]+floor(j/2))+bord,col=subcols[[i]][j+1],border=TRUE)
				text(0.3+bord+0.01,1-y.unit*(i+n.add[i]+floor(j/2))+(y.unit-bord)/2,sublabels[[i]][j+1],cex=cex.sub,adj=c(0,0.5),font=2)
			}
			
		}
	}else{
		z <- seq(sublabels[[i]][1], sublabels[[i]][2], length = length(subcols[[i]]))
		z = cbind(z,z)
		for (j in seq(1,n1[i],2)){
		image(x=seq(0.3-2*bord,0.3+2*bord,bord*4/74),y=c(1-y.unit*(i+n.add[i]+floor(j/2)),1-y.unit*(i+n.add[i]+floor(j/2))+bord/2),z,col = subcols[[i]], axes=FALSE,add=TRUE)
		text(0.3-2*bord-0.05,1-y.unit*(i+n.add[i]+floor(j/2))+(y.unit-bord)/2,sublabels[[i]][1],cex=cex.sub,adj=c(0,0.5),font=2)
		text(0.3+2*bord+0.05,1-y.unit*(i+n.add[i]+floor(j/2))+(y.unit-bord)/2,sublabels[[i]][2],cex=cex.sub,adj=c(0,0.5),font=2)
		
	}
}

}
box.coord<-list(x=c(-0.05,-0.05,0.65,0.65,-0.05),y=c(0,1.05,1.05,0,0))
lines(box.coord,xpd=TRUE)
}


############################################################################
# HISTORY:
# Jun 7, 2011
# o Created. by hui
############################################################################
