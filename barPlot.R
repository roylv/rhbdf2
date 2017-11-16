setwd("/scratch200/roeilevi/rhbdf")
#con408=read.delim("scan_ConnD14_408.txt")
#tbth=read.delim("scan_TbTh601.txt")
#bv1192=read.delim("scan_1192.txt")
#loc=matrix(c("chr7", 185,"scan_TbTh601.txt","chr2", 548, "scan_1192.txt","chr8",154, "scan_TbTh601.txt" ,"chr14",286,"scan_ConnD14_408.txt"),nrow=4,ncol=3, byrow=T) 
#loc=matrix(c("chr2", 548, "scan_1192.txt","chr2", 504, "scan_1192.txt","chr8",154, "scan_TbTh601.txt","chr7", 185,"scan_TbTh601.txt","chr7", 216,"scan_TbTh601.txt","chr14",286,"scan_ConnD14_408.txt"),nrow=4,ncol=3, byrow=T)

#df1="/Net/mus/data/www/preCC/BONE/LEVY/redHawk/MCRemoved/Trab/imputed15.txt"

#loc=matrix(c("chr7", 185,df1,"chr2", 548, df1,"chr8",154, df1 ,"chr14",286,df1),nrow=4,ncol=3, byrow=T)
source("happyEx19.preCC.R")
condensed.db = load.condensed.database("../Mar2015.MEGA+MDA+MUGA/CONDENSED/")
load("mtx2n.RData")
mtx=mtx2
db=condensed.db
coeff.names <- c("A.J", "C57", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
#loc=data.frame(loc)
#loc$X4=c("Tb.Th3D", "BVTV","Tb.NPl", "Conn.D") 
#loc$X4=c("BVTV", "BVTV","Tb.NPl","Tb.Th3D","Tb.Th3D","Conn.D")

library("RColorBrewer", "~/R")
hmcol=colorRampPalette(brewer.pal(8,"PRGn"))(10)
idxn=c(8,9, 6, 2, 1, 7)
pdf("~/public_html/rhbdf/BarAll1Feb16.pdf")

par(mfrow=c(3,2))
n=0
j=8
k=1
l=7
for (i in idxn) {
n=n+1
#  d.1 <- db$additive$chr$chr8[156][[1]]$mat

  phen.data=read.delim(mtx[[i,1]])

  agg.phen.data = aggregate( phen.data[[mtx[[i,2]]]], by=list(CC.Line=phen.data$CC.Line), FUN="mean" )
  int = intersect( db$subjects, unique(agg.phen.data$CC.Line))
  phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]
  db.use = match( phen.data$CC.Line, db$subjects)

  w=as.numeric(mtx[[i,3]])
  chr=mtx[[i,4]]

  e=db$additive$chr[[chr]][[w]]$mat
  d=e[db.use,]	
  
  #w=loc$X2[i]
  #chr=loc$X1[i]
  fit.genetic=lm(x~d, data=phen.data)
  
  #d.1 <- db$additive$chr[[loc$X1[i]]][loc$X2[i]][[1]]$mat
  #.1 = d.1[db.use,]
#$chr8[156][[1]]$mat
  #fit=lm(x~d.1, data=phen.data)
 
  coeff.vec.1 <- c(fit.genetic$coefficients[2:8])
  
	  if (i %in% c(8)){
  barplot(c(coeff.vec.1,0), names.arg=coeff.names, col=hmcol,ylab="Strain deviation", main=paste("Trl7 (BV/TV)"), cex.names=1, las=2 ); title(main="A", adj=0)}
 if (i %in% c(9)){
  barplot(c(coeff.vec.1,0), names.arg=coeff.names, col=hmcol,ylab="Strain deviation", main=paste("Trl7 (Tb.N)"), cex.names=1, las=2 ); title(main="B", adj=0)}

	if (i %in% c(6,2)){	
  barplot(c(coeff.vec.1,0), names.arg=coeff.names, col=hmcol,ylab="Strain deviation", main=paste("Trl", j, sep=""), cex.names=1, las=2 ); if (i==6) title(main="C", adj=0); if(i==2) title(main="D", adj=0); j=j+1 }
  
  else if (i %in% c(1,7))  {barplot(c(coeff.vec.1,0), names.arg=coeff.names, col=hmcol,ylab="Strain deviation", main=paste("Crl", k, sep=""), cex.names=1, las=2 ); if (i==1) title(main="E", adj=0); if(i==7) title(main="F", adj=0); k=k+1}
}
dev.off()
