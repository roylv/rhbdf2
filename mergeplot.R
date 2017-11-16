pdf("merge4Loci3(Sep15).pdf")
library(plotrix, lib="~/R")

#load("phenAll.RData")
load("mtx2n.RData")
mtx=mtx2
q=c(0.95,0.99)
for  (i in sq=c(8, 9, 6, 2, 1, 7)){
	load(paste(mtx[i,2],mtx[i,3],"N.RData", sep=""))	
#	assign(paste(mtx[[i,2]],i),read.delim(Sys.glob(file.path("merge2", paste(mtx[[i,2]],".merge",sep=""), paste(mtx[[i,2]],".chr",mtx[[i,4]],"*",sep=""))),sep=""))
        assign(paste(mtx[[i,2]],i, sep=""),read.delim(Sys.glob(file.path(paste(mtx[[i,2]],mtx[[i,3]],".merge",sep=""), paste(mtx[[i,2]],".",mtx[[i,4]],"*",sep=""))),sep=""))

	assign(paste("sp",i,sep=""),strsplit(as.character(get(paste(mtx[[i,2]],i))$sdp),".",fixed=T))
	sp_ob=paste("sp",i,sep="")
	assign(sp_ob,sapply(get(sp_ob), function(x) as.numeric(x)))
	assign(paste(sp_ob,"idx",sep=""),which(sapply(1:ncol(get(sp_ob)), function(x) sum(get(sp_ob)[,x] %in% 3:10)>0)))

	assign(paste(mtx[[i,2]],i,sep=""), `[[<-`(get(paste(mtx[[i,2]],i)), 'multiallelic', value = 0))

	ma=get(paste(sp_ob,"idx",sep=""))
	ba=get(paste(mtx[[i,2]],i,sep=""))
	ba$multiallelic[ma]=1
	assign(paste(mtx[[i,2]],i,sep=""),ba)
	load(paste("scans",mtx[[i,1]],".RData"))
	assign(paste("sc",mtx[[i,2]],sep=""), scan.data$scans[[mtx[[i,2]]]])
	gd=get(paste("sc",mtx[[i,2]],sep=""))
	thr=gd$genomewide.distribution[ floor( length(gd$genomewide.distribution)*q) ][2]
        assign(paste("q",i,sep=""),thr)
 	
	df=read.delim(paste(mtx[[i,2]],".scan.txt",sep=""))
	
	dfc=df[df$chromosome==paste("chr",mtx[[i,4]],sep=""),]
	w=which(df$from.marker==dfc$from.marker[which.max(dfc$logP)])

	assign(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""), which(get(paste(mtx[[i,2]],i,sep=""))$bp>df$from.bp[w-5] & get(paste(mtx[[i,2]],i,sep=""))$bp<df$to.bp[w+5]))

	thresh=4.9
	assign(paste("ind",mtx[[i,2]],mtx[[i,4]],sep=""), which(get(paste(mtx[[i,2]],i,sep=""))$logP.merge[get(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""))]>thresh))

	cat(i,"\n")

}

 mouseRec=read.csv("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389972/bin/supp_112.141036_TableS1.csv")
 mouseRec[,2]=mouseRec[,2]*1000

 colHs="DarkSlateGray"
 col=1#"#660033"
 colMrg="lightgray"
 colMul="IndianRed"#"#CC0000"
 par(mfrow=c(3,2))#)
 colRec="PowderBlue" #"#9999FF"
 cexMain=1.3
 lwd=0.3
 pch=19


pdf("sampMrg.pdf")

par(mai=c(.3,.15,.5,.40))
par(omi=c(.5,.4,.5,.5))


par(mfrow=c(3,2))
for (i in sq=c(8, 9, 6, 2, 1, 7)){
	mrg = get(paste(mtx[[i,2]],i,sep=""))
	int=get(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""))
	ind=get(paste("ind",mtx[[i,2]],mtx[[i,4]],sep=""))
	spidx=get(paste("sp",i,"idx",sep=""))
	Q=get(paste("q",i,sep=""))
	
	rr=mouseRec[mouseRec[,1]==paste("chr",mtx[[i,4]],sep=""),]
	min.pos=min(mrg$bp[int][ind])
	max.pos=max(mrg$bp[int][ind])
	
	keep.recomb <- subset(rr, rr[,2] > min.pos & rr[,2] < max.pos)
        fac=40	
	if (i==4) fac=5
	if (i==5) fac=3
	ylab=""
	if (i==3 | i==9) ylab="logP"
	plot(main=substitute(italic("Xxx")), cex.main=cexMain, mrg$bp[int][ind]/1e6, mrg$logP.merge[int][ind], ylim=c(0, 1.2+max( mrg$logP.merge[int][ind])), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab=ylab) #xlim=c(129.7,130.8)
#	points(mrg$bp[spidx][ind]/1e6, mrg$logP.merge[spidx][ind], pch=18, col=colMul)

	points(mrg$bp[spidx][ind][which(mrg$logP.merge[spidx][ind]>thresh)]/1e6, mrg$logP.merge[spidx][ind][which(mrg$logP.merge[spidx][ind]>thresh)],col=colMul, pch=18)
	lines(c(min(mrg$bp[int][ind]/1e6)-1,max(mrg$bp[int][ind]/1e6)+1), c(Q,Q),lty="dotted",lwd=1.3, col="black")
	lines(mrg$bp[int][ind]/1e6, mrg$logP.interval[int][ind], col=col, lwd=1)
	lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*fac, type="l", col=colRec)
	if (i!=3 & i!=4 & i!=5) axis(4, at=seq(0,12,by=2), labels=c("0",".025",".05",".75",".1",".125",".15"),las=3, cex=0.75)
	if (i==4) axis(4, at=seq(0,15,by=7.5), labels=c("0","1","2"),las=3)
	if (i==5) axis(4, at=seq(0,12,by=6), labels=c("0","2","4"),las=3)
	if (i==4 | i==10) mtext(text=expression("4N"[e]*"r/Kb"), side=4, line=3, cex=0.75)
        if (i==3 | i==9) mtext(text="logP", side=2, line=3, cex=0.75)



}
dev.off()





#tnmrg=read.delim("Tb.NPl.merge/Tb.NPl.chr8:4500000-44400000.merge.txt", sep=" ")
#thmrg=read.delim("Tb.Th3D.merge/Tb.Th3D.chr7:62200000-75300000.merge.txt", sep=" ")
#bmrg=read.delim("BVTV.merge/BVTV.chr2:1.22e+08-130700000.merge.txt", sep=" ") 
#cmrg=read.delim("Conn.D.merge/Conn.D.chr14:69300000-91500000.merge.txt", sep=" ")
#spn=strsplit(as.character(tnmrg$sdp), ".", fixed=T)
#spn=sapply(spn, function(x) as.numeric(x));spnidx=which(sapply(1:ncol(spn), function(x) sum(spn[,x] %in% 3:10)>0));tnmrg$multiallelic=0; tnmrg$multiallelic[spnidx]=1
#sph=strsplit(as.character(thmrg$sdp), ".", fixed=T)
#sph=sapply(sph, function(x) as.numeric(x));sphidx=which(sapply(1:ncol(sph), function(x) sum(sph[,x] %in% 3:10)>0)); thmrg$multiallelic=0; thmrg$multiallelic[sphidx]=1
#spc=strsplit(as.character(cmrg$sdp), ".", fixed=T)
#spc=sapply(spc, function(x) as.numeric(x));spcidx=which(sapply(1:ncol(spc), function(x) sum(spc[,x] %in% 3:10)>0)); cmrg$multiallelic=0; cmrg$multiallelic[spcidx]=1
#spb=strsplit(as.character(bmrg$sdp), ".", fixed=T)
#spb=sapply(spb, function(x) as.numeric(x));spbidx=which(sapply(1:ncol(spb), function(x) sum(spb[,x] %in% 3:10)>0)); bmrg$multiallelic=0; bmrg$multiallelic[spbidx]=1

#*******LOAD SCANS FOR THRESHOLD*****#
#q=c(0.95,0.99)
#load("TbTh601.RData")
#scTh=scan.data$scans$Tb.Th3D
#scTn=scan.data$scans$Tb.NPl
#load("scans_1192.RData")
#scTb=scan.data$scans$BVTV
#load("ConnD14_408.RData")
#scCon=scan.data$scans$Conn.D

#qH= scTh$genomewide.distribution[ floor( length(scTh$genomewide.distribution)*q) ][2]
#qN= scTn$genomewide.distribution[ floor( length(scTn$genomewide.distribution)*q) ][2]
#qB= scTb$genomewide.distribution[ floor( length(scTb$genomewide.distribution)*q) ][2]
#qC= scCon$genomewide.distribution[ floor( length(scCon$genomewide.distribution)*q) ][2]






#wn=5123
#w=1452
#wc=8602
#wh=4535
#wh2=4566
#thresh=4.5
#lwd=0.3
#pch=19



#intb=which(bmrg$bp>bv$from.bp[w-5] & bmrg$bp<bv$to.bp[w+1])
#indb=which(bmrg$logP.merge[intb]>3.5)
#intc=which(cmrg$bp>con$from.bp[wc-9] & cmrg$bp<con$to.bp[wc+15])
#indc=which(cmrg$logP.merge[intc]>thresh)

#inth=which(thmrg$bp>th$from.bp[wh+1] & thmrg$bp<th$to.bp[wh+6])
#indh=which(thmrg$logP.merge[inth]>thresh)


#inth2=which(thmrg$bp>th$from.bp[wh2-3] & thmrg$bp<th$to.bp[wh2+3])
#indh2=which(thmrg$logP.merge[inth2]>thresh)

#intn=which(tnmrg$bp>tn$from.bp[wn-10] & tnmrg$bp<tn$to.bp[wn+6])
#indn=which(tnmrg$logP.merge[intn]>thresh)
#intb2=which(bmrg$bp>bv$from.bp[w-44-5] & bmrg$bp<bv$to.bp[w-44+5])
#indb2=which(bmrg$logP.merge[intb2]>thresh)

cn2=url("http://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/genetic_map_chr2.txt")
rec2=read.delim(cn2, sep="")
cn14=url("http://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/genetic_map_chr14.txt")
rec14=read.delim(cn14, sep="")
cn8=url("http://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/genetic_map_chr8.txt")
rec8=read.delim(cn8, sep="")
cn7=url("http://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/genetic_map_chr7.txt")
rec7=read.delim(cn7, sep="")

mouseRec=read.csv("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389972/bin/supp_112.141036_TableS1.csv")
mouseRec[,2]=mouseRec[,2]*1000
rr2=mouseRec[mouseRec[,1]=="chr2",]
rr7=mouseRec[mouseRec[,1]=="chr7",]
rr8=mouseRec[mouseRec[,1]=="chr8",]
rr14=mouseRec[mouseRec[,1]=="chr14",]
hs=read.csv("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389972/bin/supp_112.141036_TableS2.csv")
hs2=hs[hs[,1]==2,]
hs7=hs[hs[,1]==7,]
hs8=hs[hs[,1]==8,]
hs14=hs[hs[,1]==14,]


 colHs="DarkSlateGray"
 col=1#"#660033"
 colMrg="lightgray"
 colMul="IndianRed"#"#CC0000"
 par(mfrow=c(3,2))#)
 colRec="PowderBlue" #"#9999FF"
 par(mai=c(.3,.15,.5,.40))
 par(omi=c(.5,.4,.5,.5))
 cexMain=1.3
 plot(main=substitute(italic("Trl1")), cex.main=cexMain, bmrg$bp[intb][indb]/1e6, bmrg$logP.merge[intb][indb], ylim=c(0,8.5), lwd=lwd, col=colMrg, xlim=c(129.7,130.8), pch=pch, xlab="Mb", ylab="logP")


#**********RECOMBINATION RATES*************#
min.pos=129800000
max.pos=130700000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

keep.recomb <- subset(rr2, rr2[,2] > min.pos & rr2[,2] < max.pos)
lines(keep.recomb[,2]/1e6,  keep.recomb[,3]* 64, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)
keep.recombhs <- subset(hs2, hs2[,2] > min.pos & hs2[,2] < max.pos)
#points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=5,col=colHs,bg=224,cex=0.7)

lines(c(min(bmrg$bp[intb][indb]/1e6)-1,max(bmrg$bp[intb][indb]/1e6)+1), c(qB,qB),lty="dotted",lwd=1.3, col="black")



axis(4, at=seq(0,8,by=4), labels=c("0","0.05","0.1"),las=3) 
#**********RECOMBINATION RATES*************#

ind=which(bmrg$logP.merge[spbidx]>thresh & bmrg$bp[spbidx]>129100000)
 points(bmrg$bp[spbidx][ind]/1e6, bmrg$logP.merge[spbidx][ind], pch=18, col=colMul)

 lines(bmrg$bp[intb][indb]/1e6, bmrg$logP.interval[intb][indb], ylim=c(0,8.5), col=col, lwd=1)
        rect(130.576459,0.5,130.587523,1.1, col=3)
        rect(130.571252,1.1,130.582053,1.7, col=1)
#	rect(130.27,0.5,130.28,1.1,col=2)
	legend("topleft", legend=c("MA", "BA","Oxt","Avp"), col=c(colMul,colMrg,3,1), pch=c(18,1,15,15),cex=0.71)
	title(main="A", adj=0)


 plot(main=substitute(italic("Trl2")), cex.main=cexMain,bmrg$bp[intb2][indb2]/1e6, bmrg$logP.merge[intb2][indb2], ylim=c(0,8.5), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(121.9,max(bmrg$bp[intb2][indb2]/1e6)))


#**********RECOMBINATION RATES*************#
min.pos=122000000
max.pos=122600000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#keep.recomb <- subset(rec2, rec2[,1] > min.pos & rec2[,1] < max.pos)
#lines(keep.recomb[,1]/1e6,  keep.recomb[,2] / 10, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)

 keep.recomb <- subset(rr2, rr2[,2] > min.pos & rr2[,2] < max.pos)
 lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*640, type="l", col=colRec)

lines(c(min(bmrg$bp[intb2][indb2]/1e6)-1,max(bmrg$bp[intb2][indb2]/1e6)+1), c(qB,qB),lty="dotted",lwd=1.3, col="black")

keep.recombhs <- subset(hs2, hs2[,2] > min.pos & hs2[,2] < max.pos)
#points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=5,col="maroon",cex=0.7,bg="blue")



#axis(4, at=seq(0,8,by=2), labels=c("0","4","8","12","18"),las=1)
axis(4, at=seq(0,8,by=4), labels=c("0",expression("1e-5"),"1e-4"),las=3)

#**********RECOMBINATION RATES*************#



 ind=which(bmrg$logP.merge[spbidx]>thresh & bmrg$bp[spbidx]>121.9*1e6)
 points(bmrg$bp[spbidx][ind]/1e6, bmrg$logP.merge[spbidx][ind], pch=18, col=colMul)

 lines(bmrg$bp[intb2][indb2]/1e6, bmrg$logP.interval[intb2][indb2], ylim=c(0,8.5), col=col, lwd=1)
         rect(122142790/1e6,0.5,122155639/1e6,1.1,col="purple")
        legend("topleft", legend=c("MA", "BA", "B2m"), col=c(colMul,colMrg,"purple"), pch=c(18,1,15), cex=0.71)
        title(main="B", adj=0)





 plot(main=substitute(bold(italic("Trl3"))), cex.main=cexMain,tnmrg$bp[intn][indn]/1e6, tnmrg$logP.merge[intn][indn], ylim=c(0,8.5), lwd=lwd, col=colMrg, xlim=c(36.1,max(tnmrg$bp[intn][indn]/1e6)), pch=pch, xlab="Mb",ylab="logP")

#**********RECOMBINATION RATES*************#
min.pos=36700000
max.pos=40730000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#keep.recomb <- subset(rec8, rec8[,1] > min.pos & rec8[,1] < max.pos)
#lines(keep.recomb[,1]/1e6,  keep.recomb[,2] / 10, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)
keep.recomb <- subset(rr8, rr8[,2] > min.pos & rr8[,2] < max.pos)
lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*16, type="l", col=colRec)

keep.recombhs <- subset(hs8, hs8[,2] > min.pos & hs8[,2] < max.pos)
#points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=22,bg="blue",cex=0.5)



lines(c(min(tnmrg$bp[intn][indn]/1e6),max(tnmrg$bp[intn][indn]/1e6)+1), c(qN,qN),lty="dotted",lwd=1.3, col="black")




#axis(4, at=seq(0,8,by=2), labels=c("0","20","40","60","80"),las=1)
axis(4, at=seq(0,8,by=4), labels=c("0","0.25","0.5"),las=3)

mtext("logP", side=2, line=3, cex=0.75)

#**********RECOMBINATION RATES*************#







 ind=which(tnmrg$logP.merge[spnidx]>thresh & tnmrg$bp[spnidx]>36.8*1e6)
 points(tnmrg$bp[spnidx][ind]/1e6, tnmrg$logP.merge[spnidx][ind], pch=18, col=colMul)
 lines(tnmrg$bp[intn][indn]/1e6, tnmrg$logP.interval[intn][indn], ylim=c(0,8.5), col=col, lwd=1)
 #ablineclip(h=4.9, x1=min(tnmrg$bp[intn][indn]/1e6),x2=max(tnmrg$bp[intn][indn]/1e6), lty=2, col="gray", lwd=1.7)
	rect(37.517382,0.5,38.666439,1.1,col="orange")
	rect(40.274167,0.5,40.313316,1.1,col="coral")
	rect(40.487547,0.5,40.520789,1.1,col="8")
	legend("topleft", legend=c("MA", "BA", "Sgcz","Fgf20","Cnot7"), col=c(colMul, colMrg, "orange","coral","8"), pch=c(18,1,15,15,15), cex=0.71)
	title(main="C", adj=0)





 plot(main=substitute(italic("Trl4")), cex.main=cexMain, thmrg$bp[inth][indh]/1e6, thmrg$logP.merge[inth][indh], ylim=c(0,9.5), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(64.5, max(thmrg$bp[inth][indh]/1e6)))


#**********RECOMBINATION RATES*************#
min.pos=64500000
max.pos=66400000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#keep.recomb <- subset(rec7, rec7[,1] > min.pos & rec7[,1] < max.pos)
#lines(keep.recomb[,1]/1e6,  keep.recomb[,2] / 10, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)

keep.recomb <- subset(rr7, rr7[,2] > min.pos & rr7[,2] < max.pos)
lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*16, type="l", col=colRec)
lines(c(min(thmrg$bp[inth][indh]/1e6)-1,max(thmrg$bp[inth][indh]/1e6)+1), c(qH,qH),lty="dotted",lwd=1.4, col="black")

keep.recombhs <- subset(hs7, hs7[,2] > min.pos & hs7[,2] < max.pos)
#points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=22,bg="blue",cex=0.5)



#axis(4, at=seq(0,8,by=2), labels=c("0","20","40","60","80"),las=1)
axis(4, at=seq(0,8,by=4), labels=c("0","0.25","0.5"),las=3)

mtext(text=expression("4N"[e]*"r/Kb"), side=4, line=3, cex=0.75) 
#mtext("logP", side=2, line=1, cex=0.75)

#**********RECOMBINATION RATES*************#





 ind=which(thmrg$logP.merge[sphidx]>thresh & thmrg$bp[sphidx]>64.9*1e6)
 points(thmrg$bp[sphidx][ind]/1e6, thmrg$logP.merge[sphidx][ind], pch=18, col=colMul)

 lines(thmrg$bp[inth][indh]/1e6, thmrg$logP.interval[inth][indh], ylim=c(0,9.5), col=col, lwd=1)

# 	 rect(64346758/1e6,1.1,64374095/1e6,1.7,col=101)#fan1
#	 rect(64376576/1e6,0.5,64392268/1e6,1.1,col=106) #mphosph10
	 
#	 rect(65856837/1e6,0.5,66055280/1e6,1.1,col="blue")
#	 rect(64501706/1e6,0.5,64753870/1e6,1.1,col="brown")#Apba2
#	 rect(64346758/1e6,1.1,64374095/1e6,1.7,col=101)#fan1
#	 rect(64392645/1e6,1.1,64412121/1e6,1.7,col=102)#mcee
#	 rect(64376576/1e6,0.5,64392268/1e6,1.1,col=106) #mphosph10

#	 rect(64867052/1e6,0.5,64872997/1e6,1.1,col=107)#Ndnl2
	 rect(65856837/1e6,0.5,66055280/1e6,1.1,col="Salmon")
	# rect(65296165/1e6,0.5,65371239/1e6,1.1,col=570)#Tjp1


	
#	 legend("topleft", legend=c("MA", "BA", "Fan1", "Mphosph10","Apba2", "Ndnl2", "Pcsk6"), col=c(colMul, colMrg, 101, 106, "brown", 107, "blue"), pch=c(18,1,15,15,15,15,15,15,15), cex=0.71)

 legend("topleft", legend=c("MA", "BA", "Pcsk6"), col=c(colMul, colMrg, "Salmon"), pch=c(18,1,15), cex=0.71)

	 title(main="D", adj=0)
#	 abline(v=c(64.5,64.75,64.86,64.87, 64.2, 64.3, 64.34, 64.37))
#	 abline(v=c(64.39, 64.41), col="red")



#**************#
 plot(main=substitute(italic("Trl5")), cex.main=cexMain, thmrg$bp[inth2][indh2]/1e6, thmrg$logP.merge[inth2][indh2], ylim=c(0,9.5), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(72.9, max(thmrg$bp[inth2][indh2]/1e6)))

#**********RECOMBINATION RATES*************#
min.pos=73100000
max.pos=74000000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#keep.recomb <- subset(rec7, rec7[,1] > min.pos & rec7[,1] < max.pos)
#lines(keep.recomb[,1]/1e6,  keep.recomb[,2] / 10, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)

 keep.recomb <- subset(rr7, rr7[,2] > min.pos & rr7[,2] < max.pos)
 lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*64, type="l", col=colRec)
lines(c(min(thmrg$bp[inth2][indh2]/1e6)-1,max(thmrg$bp[inth2][indh2]/1e6)+1), c(qH,qH),lty="dotted",lwd=1.3, col="black")

keep.recombhs <- subset(hs7, hs7[,2] > min.pos & hs7[,2] < max.pos)
#points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=22,bg="blue",cex=0.5)



#axis(4, at=seq(0,8,by=2), labels=c("0","20","40","60","80"),las=1)
axis(4, at=seq(0,8,by=4), labels=c("0","0.05","0.1"),las=3)

#**********RECOMBINATION RATES*************#



 ind=which(thmrg$logP.merge[sphidx]>4 & thmrg$bp[sphidx]>73.1*1e6)
 points(thmrg$bp[sphidx][ind]/1e6, thmrg$logP.merge[sphidx][ind], pch=18, col=colMul)

 lines(thmrg$bp[inth2][indh2]/1e6, thmrg$logP.interval[inth2][indh2], ylim=c(0,9.5), col=col, lwd=1)



        rect(73375520/1e6,0.5,73419899/1e6,1.1, col=14)#rgma
        rect(73426633/1e6,1.1,73541746/1e6,1.7, col=15)#chd2

	legend("topleft", legend=c("MA", "BA","Rgma", "Chd2"), col=c(colMul, colMrg, 14,15), pch=c(18,1,15,15), cex=0.71)

	title(main="E", adj=0)


#*************#






 plot(main=substitute(italic("Trl6")), cex.main=cexMain,cmrg$bp[intc][indc]/1e6, cmrg$logP.merge[intc][indc], ylim=c(0,8.5),lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab="logP", xlim=c(72.42,max(cmrg$bp[intc][indc]/1e6)))


 ind=which(cmrg$logP.merge[spcidx]>thresh & cmrg$bp[spcidx]>73000000)
 min.pos=73000000
 max.pos=77000000
#*******RECOMBINATION RATES******#
 size.pos <- max.pos - min.pos
 center.pos <- min.pos + ( size.pos / 2 )
 center.100kb.pos <- round(center.pos / 100000) * 100000
 offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

 #keep.recomb <- subset(rec14, rec14[,1] > min.pos & rec14[,1] < max.pos) 
 #lines(keep.recomb[,1]/1e6,  keep.recomb[,2] / 10, type="l", col=colRec)#, lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset

 keep.recomb <- subset(rr14, rr14[,2] > min.pos & rr14[,2] < max.pos)
 lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*3.7, type="l", col=colRec)
 
 keep.recombhs <- subset(hs14, hs14[,2] > min.pos & hs14[,2] < max.pos)
 #points(keep.recombhs[,2]/1e6,rep(2.5,nrow(keep.recombhs)),pch=22,bg="blue",cex=0.5)
 #points(keep.recombhs[,3]/1e6,rep(2.5,nrow(keep.recombhs)),pch=22,bg="red",cex=0.5)


 axis(4, at=seq(0,8,by=4), labels=c("0","1","2"),las=1)

 
#*******RECOMBINATION RATES******#

 lines(c(min(cmrg$bp[intc][indc]/1e6),max(cmrg$bp[intc][indc]/1e6)+1), c(qC,qC),lty="dotted",lwd=1.3, col="black")
 

# ind2=which(cmrg$bp>72420000)
# ind3=which(ind %in% ind2)
 points(cmrg$bp[spcidx][ind]/1e6, cmrg$logP.merge[spcidx][ind], pch=18, col=colMul)
 lines(cmrg$bp[intc][indc]/1e6, cmrg$logP.interval[intc][indc], ylim=c(0,8.5), col=col, lwd=1)
#        rect(75237306/1e6,0.5,75288508/1e6,1.1, col=3)#cpb2 
#        rect(73194749/1e6,0.5,73330793/1e6,1.1, col=4)#rb1
#        rect(74754871/1e6,0.5,74952722/1e6,1.1, col=5)#lrch1
#        rect(74727348/1e6,1.1,74754650/1e6,1.7, col=6)#esd

        rect(73194749/1e6,0.5,73330793/1e6,1.1, col=4)#rb1
        rect(74727348/1e6,1.1,74754650/1e6,1.7, col=6)#esd
        rect(74754871/1e6,0.5,74952722/1e6,1.1, col=5)#lrch1
        rect(75237306/1e6,0.5,75288508/1e6,1.1, col=3)#cpb2



        legend("topleft", legend=c("MA", "BA","Rb1", "Esd", "Lrch1", "Cpb2"), col=c(colMul, colMrg, 4,6,5,3), pch=c(18,1,15,15,15,15), cex=0.71)
        title(main="F", adj=0)
	

 mtext("Mb", side=1, line=3, adj=-.145, cex=.75)

 dev.off()
