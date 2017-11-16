library(plotrix)
getwd()
#load("phenAll.RData")
load("mtx2n.RData")
mtx=mtx3
q=c(0.95,0.99)
sq=c(8, 9, 6, 2, 1, 7)
#sq=c(8,6)
for  (i in sq){
        load(paste(mtx[i,2],mtx[i,3],"N.RData", sep=""))
	
#       assign(paste(mtx[[i,2]],i),read.delim(Sys.glob(file.path("merge2", paste(mtx[[i,2]],".merge",sep=""), paste(mtx[[i,2]],".chr",mtx[[i,4]],"*",sep=""))),sep=""))
        assign(paste(mtx[[i,2]],i, sep=""),read.delim(Sys.glob(file.path(paste(mtx[[i,2]],mtx[[i,3]],".merge",sep=""), paste(mtx[[i,2]],".",mtx[[i,4]],"*",sep=""))),sep=""))
	
        assign(paste("sp",i,sep=""),strsplit(as.character(get(paste(mtx[[i,2]],i,sep=""))$sdp),".",fixed=T))
	
        sp_ob=paste("sp",i,sep="")
        assign(sp_ob,sapply(get(sp_ob), function(x) as.numeric(x)))
        assign(paste(sp_ob,"idx",sep=""),which(sapply(1:ncol(get(sp_ob)), function(x) sum(get(sp_ob)[,x] %in% 3:10)>0)))

        assign(paste(mtx[[i,2]],i,sep=""), `[[<-`(get(paste(mtx[[i,2]],i, sep="")), 'multiallelic', value = 0))

        ma=get(paste(sp_ob,"idx",sep=""))
        ba=get(paste(mtx[[i,2]],i,sep=""))
        ba$multiallelic[ma]=1
        assign(paste(mtx[[i,2]],i,sep=""),ba)
        #load(paste("scans",mtx[[i,1]],".RData"))
        assign(paste("sc",mtx[[i,2]],sep=""), scan.data$scans[[mtx[[i,2]]]])
        gd=get(paste("sc",mtx[[i,2]],sep=""))
        thr=gd$genomewide.distribution[ floor( length(gd$genomewide.distribution)*q) ][2]
        assign(paste("q",i,sep=""),thr)

        #df=read.delim(paste(mtx[[i,2]],".scan.txt",sep=""))
	df=scan.data$scans[[mtx[i,2]]]$result
        dfc=df[df$chromosome==paste(mtx[[i,4]],sep=""),]
        w=which(df$from.marker==dfc$from.marker[which.max(dfc$logP)])

        assign(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""), which(get(paste(mtx[[i,2]],i,sep=""))$bp>df$from.bp[w-5] & get(paste(mtx[[i,2]],i,sep=""))$bp<df$to.bp[w+5]))

        thresh=5
        assign(paste("ind",mtx[[i,2]],mtx[[i,4]],sep=""), which(get(paste(mtx[[i,2]],i,sep=""))$logP.merge[get(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""))]>thresh))

        cat(i,"\n")

}




# mouseRec=read.csv("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389972/bin/supp_112.141036_TableS1.csv")
mouseRec=read.csv("recomb.csv") 
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


pdf("~/public_html/rhbdf/mrgAug16.2.pdf")

par(mai=c(.3,.3,.5,.40))
par(omi=c(.5,.3,.5,.5))
par(oma=c(4,2,1,1))
par(xpd=T)

Xxx=c("Trl7 (BV/TV)","Trl7 (Tb.N)", "Trl8", "Trl9", "Crl1", "Crl2")

cxl=0.9
j=1
#par(xpd=T)
par(mfrow=c(3,2))
for (i in sq){
        mrg = get(paste(mtx[[i,2]],i,sep=""))

        int=get(paste("int",mtx[[i,2]],mtx[[i,4]],sep=""))
        ind=get(paste("ind",mtx[[i,2]],mtx[[i,4]],sep=""))
        spidx=get(paste("sp",i,"idx",sep=""))
        Q=get(paste("q",i,sep=""))

        rr=mouseRec[mouseRec[,1]==paste(mtx[[i,4]],sep=""),]
        min.pos=min(mrg$bp[int][ind])
        max.pos=max(mrg$bp[int][ind])

        keep.recomb <- subset(rr, rr[,2] > min.pos & rr[,2] < max.pos)
        if (i %in% c(8,9)) {fac=40
	}else fac=6
        if (i==4) fac=5
        if (i==5) fac=3
        ylab=""
        #if (i==3 | i==9) ylab="logP"
	if (i==8) {strt=0.3; endd=0.2}
	if (i==9) {strt=0.3; endd=0.9}
	if (i==6) {strt=1.25; endd=-0.03
	
	}else {endd=0; strt=0}

	mrg2=mrg[int,][ind,][which(mrg$bp[int][ind]/1e6>((mrg$bp[int][ind]/1e6)[1]+strt) & mrg$bp[int][ind]/1e6<(tail(mrg$bp[int][ind]/1e6, n=1)-endd)),]
        min.pos=min(mrg2$bp)
        max.pos=max(mrg2$bp)
	keep.recomb <- subset(rr, rr[,2] > min.pos & rr[,2] < max.pos)



#	plot(main=Xxx[j], cex.main=cexMain, mrg$bp[int][ind]/1e6, mrg$logP.merge[int][ind], ylim=c(0, 1.2+max( mrg$logP.merge[int][ind])), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab=ylab, xlim=c((mrg$bp[int][ind]/1e6)[1]+strt,tail(mrg$bp[int][ind]/1e6, n=1)-endd))
	plot(main=Xxx[j], cex.main=cexMain, mrg2$bp/1e6, mrg2$logP.merge, ylim=c(0, 1.2+max( mrg2$logP.merge)), lwd=lwd, col=colMrg, pch=pch, xlab="Mb", ylab=ylab,frame=F) #xlim=c((mrg$bp[int][ind]/1e6)[1]+strt,tail(mrg$bp[int][ind]/1e6, n=1)-endd))
	 

	j=j+1
#       points(mrg$bp[spidx][ind]/1e6, mrg$logP.merge[spidx][ind], pch=18, col=colMul)

        points(mrg2$bp[spidx][which(mrg2$logP.merge[spidx]>thresh)]/1e6, mrg2$logP.merge[spidx][which(mrg2$logP.merge[spidx]>thresh)],col=colMul, pch=18)

#	 lines(c(min(mrg$bp[int][ind]/1e6)-1,max(mrg$bp[int][ind]/1e6)+1), c(Q,Q),lty="dotted",lwd=1.3, col="black")

        lines(c(min(mrg2$bp/1e6),max(mrg2$bp/1e6)), c(Q,Q),lty="dotted",lwd=1.3, col="black")
        lines(mrg2$bp/1e6, mrg2$logP.interval, col=col, lwd=1)
        lines(keep.recomb[,2]/1e6,  keep.recomb[,3]*fac, type="l", col=colRec, lwd=1.4)
        


	if (i==8) axis(4, at=seq(0,12,by=4), labels=round(seq(0,0.104,0.034),2),las=3, cex=0.75)
	if (i==9) axis(4, at=seq(0,10,by=2), labels=round(seq(0,0.08,0.016),2),las=3, cex=0.75)


        if (i==6 | i==1) axis(4, at=seq(0,10,by=2), labels=round(seq(0,1.825,0.365),2),las=3)
        
	if (i==2) axis(4, at=seq(0,12,by=2), labels=c(seq(0,3.5, 0.7),las=3))
        if (i==2) mtext(text=expression("4N"[e]*"r/Kb"), side=4, line=3, cex=0.75)
        if (i==6) mtext(text="logP", side=2, line=3, cex=0.75)

	if (i==7) axis(4, at=seq(0,12,by=4), labels=round(seq(0,2.34,0.78),2),las=3, cex=0.75)
	
	#legends

	if (i==8) {
        rect(116598165/1e6,0.5,116627019/1e6,1.5, col="#996633", border="#996633", bty="n")#chd2
        rect(116537740/1e6,0.5,116581447/1e6,1.5, col="orange", border="orange", bty="n")#chd2
	rect(116586432/1e6,0.5,116597680/1e6,1.5, col="orange", border="orange", bty="n")
        rect(116849901/1e6,0.5,116853094 /1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2
        rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF6666", border="#FF6666", bty="n")#chd2


        #rect(116849901/1e6,0.5,116853094 /1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2
        #rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF6666", border="#FF6666", bty="n")#chd2


#	abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#	legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")
	 legend("topleft", horiz=T, bty="n",  legend=c("Ube2o","Rhdf2"), col=c("orange", "#996633"), pch=c(15,15), cex=cxl, bg="white", box.col="white", inset=c(0,-0.15))

	title(main="A", adj=0)
	}



        if (i==9) {
	
	rect(116598165/1e6,0.5,116627019/1e6,1.5, col="#996633", border="#996633", bty="n")#chd2
        rect(116537740/1e6,0.5,116581447/1e6,1.5, col="orange", border="orange", bty="n")#chd2
        rect(116849901/1e6,0.5,116853094 /1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2
        rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF6666", border="#FF6666", bty="n")#chd2
	rect(116586432/1e6,0.5,116597680/1e6,1.5, col="orange", border="orange", bty="n")




#       abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#       legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")
         legend("topleft", horiz=T, bty="n", legend=c("Ube2o","Rhdf2", "St6galnac2","Srsf2" ), col=c("orange", "#996633", "#FF6666", "#FF0000"), pch=c(15,15, 15, 15), cex=cxl, bg="white", box.col="white", inset=c(0,-0.15))

        title(main="B", adj=0)
        }





        if (i==6) {
        rect(117757836/1e6,0.5,117765648/1e6,1.5, col="orange", border="orange", bty="n")#chd2
        rect(118136957/1e6,0.5,118180043/1e6,1.5, col="#660066", border="#660066", bty="n")#chd2
#        rect(116849901/1e6,0.5,116853094/1e6,1.5, col="#FF6666", border="#FF6666", bty="n")#chd2
#        rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2


#       abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#       legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")
         legend("topleft", horiz=T, bty="n", legend=c("Klf17", "Kdm4a" ), col=c("orange", "#660066"), pch=c(15,15), cex=cxl, bg="white", box.col="white", inset=c(0,-0.15))

        title(main="C", adj=0)
        }


	if(i==2){


        rect(106452523/1e6,0.5,106458166/1e6,1.5, col="orange", border="orange", bty="n")#chd2
        rect(106616739/1e6,0.5,106697287/1e6,1.5, col="#CC99FF", border="#CC99FF", bty="n")#chd2
#        rect(116849901/1e6,0.5,116853094/1e6,1.5, col="#FF6666", border="#FF6666", bty="n")#chd2
#        rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2


#       abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#       legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")
         legend("topleft", horiz=T, bty="n", legend=c("Barhl2", "Zfp644" ), col=c("orange", "#CC99FF"), pch=c(15,15), cex=cxl, bg="white", box.col="white", inset=c(0,-0.15))
        title(main="D", adj=0)
	}





        if (i==1) {
        rect(9269293/1e6,0.5,9451691/1e6,1.5, col="orange", border="orange", bty="n")#chd2
        rect(9844372/1e6,0.5,9862345/1e6,1.5, col="#996633", border="#996633", bty="n")#chd2
        rect(9448069/1e6,1.5,9669344/1e6,2.5, col="#FF6666", border="#FF6666", bty="n")#chd2
#        rect(116675303/1e6,0.5,116694868/1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2


#       abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#       legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")


         legend("topleft", horiz=T, bty="n", legend=c("Clvs1", "Asph", "Gdf6" ), col=c("orange", "#FF6666", "#996633"), pch=c(15,15, 15), cex=cxl, bg="white", box.col="white", inset=c(0,-0.15))

        title(main="E", adj=0)
        }


        if (i==7) {
#        rect(94612757/1e6,0.5,94658872/1e6,1.5, col="orange", border="orange", bty="n")#chd2
        rect(95323525/1e6,0.5,95357202/1e6,1.5, col="#660066", border="#660066", bty="n")#chd2
        rect(95499256/1e6,1.5,95509362/1e6,2.5, col="#FF6666", border="#FF6666", bty="n")#chd2
        rect(96525172/1e6,0.5,96529210/1e6,1.5, col="#FF0000", border="#FF0000", bty="n")#chd2
        rect(97203662/1e6,0.5,97297917/1e6,1.5, col="#9999FF", border="#9999FF", bty="n")#chd2
	rect(97158777/1e6,0.5,97177299/1e6,1.5, col="#6699CC", border="#6699CC", bty="n")#chd2
        rect(98013527/1e6,0.5,98150361/1e6,1.5, col="#CC99FF", border="#CC99FF", bty="n")#chd2
#	rect(98013527/1e6,0.5,98150361/1e6,1.5, col="#CC99FF", border="#CC99FF", bty="n")#chd2

#       abline(v=116598165/1e6, lwd=0.1, color="lightgray"); abline(v=116627019/1e6, lwd=0.1, color="lightgray")
#       legend("topright", horiz=T, bty="n", legend=c("MA", "BA","Rhdf2"), col=c(colMul, colMrg, "orange"), pch=c(18,1,15), cex=0.87, bg="white")
         legend("topleft", horiz=T, bty="n", legend=c("Hfe2","Acp6","Bcl9", "Notch2" ), col=c("#FF0000","#6699CC","#9999FF", "#CC99FF"), pch=c(15,15, 15, 15), cex=cxl, bg="white", box.col="white",inset=c(0,-0.15))

        title(main="F", adj=0)
        }





}
mtext("Mb", side=1, line=3, adj=-.145, cex=.75)

dev.off()


