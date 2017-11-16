library(plotrix)

#load("phenAll.RData")
load("mtx2n.RData")
mtx=mtx2
q=c(0.95,0.99)
sq=c(8, 9, 6, 2, 1, 7)
for  (i in sq){
        load(paste(mtx[i,2],mtx[i,3],"N.RData", sep=""))
#       assign(paste(mtx[[i,2]],i),read.delim(Sys.glob(file.path("merge2", paste(mtx[[i,2]],".merge",sep=""), paste(mtx[[i,2]],".chr",mtx[[i,4]],"*",sep=""))),sep=""))
        assign(paste(mtx[[i,2]],i, sep=""),read.delim(Sys.glob(file.path(paste(mtx[[i,2]],mtx[[i,3]],".merge",sep=""), paste(mtx[[i,2]],".",mtx[[i,4]],"*",sep=""))),sep=""))

        assign(paste("sp",i,sep=""),strsplit(as.character(get(paste(mtx[[i,2]],i, sep=""))$sdp),".",fixed=T))
        sp_ob=paste("sp",i,sep="")
        assign(sp_ob,sapply(get(sp_ob), function(x) as.numeric(x)))
        assign(paste(sp_ob,"idx",sep=""),which(sapply(1:ncol(get(sp_ob)), function(x) sum(get(sp_ob)[,x] %in% 3:10)>0)))

        assign(paste(mtx[[i,2]],i,sep=""), `[[<-`(get(paste(mtx[[i,2]],i)), 'multiallelic', value = 0))

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


pdf("~/public_html/rhbdf/sampMrg.pdf")

par(mai=c(.3,.15,.5,.40))
par(omi=c(.5,.4,.5,.5))




par(mfrow=c(3,2))
for (i in 1:12){
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
#       points(mrg$bp[spidx][ind]/1e6, mrg$logP.merge[spidx][ind], pch=18, col=colMul)

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


