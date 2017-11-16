boxplot(df2[[1]]~df2[[2]], notch=T, lwd=2, pch=20, names=c("TT","T/C","CC"), ylab="BV/TV", col="grey", axes=F)
axis(side=2, at=c(0.05,0.10,0.15,0.20,0.25), labels=c(0.05,0.10,0.15,0.20,0.25))
axis(side=1, at=c(0,1,2,3,4), labels=c("","TT","T/C","CC",""))    
stripchart(df2[[1]]~df2[[2]], add=T, vertical=T, pch=20, col=c("orange","magenta","blue"), method="jitter", lwd=4)



vc2=na.omit(df[df$CC.Line %in% bvs[[1]][gc], c("BV.TV","SEX")])



 stripchart(df3[df3$Sex=="M",1]~df3[df3$Sex=="M",3], add=T, vertical=T, pch=20, col=c("blue"), method="jitter", lwd=4, jitter=0.2)
 stripchart(df3[df3$Sex=="F",1]~df3[df3$Sex=="F",3], add=T, vertical=T, pch=17, col=c("pink"), method="jitter", lwd=4, jitter=0.2)


alna=dfn[which(dfn$CC.Line %in% bvs[[1]][which(bvs[[2]]==bvs[[3]] & bvs[[2]]=="T")]), c("Tb.NPl", "Sex")]
alnc=dfn[which(dfn$CC.Line %in% bvs[[1]][which(bvs[[2]]==bvs[[3]] & bvs[[2]]=="C")]), c("Tb.NPl", "Sex")]

alnb=dfn[which(dfn$CC.Line %in% bvs[[1]][which(bvs[[2]]!=bvs[[3]])]), c("Tb.NPl", "Sex")]


 boxplot(df4[[1]]~df4$at, col=c("blue","pink"), axes=F, outline=F, notch=F)
 axis(side=1, at=c(0,1.5,3.5,5.5,7), labels=c("","TT","TC","CC",""))  
 stripchart(df4[[1]]~df4$at, col=c("black"), vertical=T, add=T, method="jitter", jitter=0.2, pch=20)


al=bvs
alPhen=function(phen, al){

 dfx=df[which(!is.na(df[[phen]])),]
 ag=aggregate( dfx[[phen]], by=list(CC.Line=dfx$CC.Line), FUN="mean")
 al=al[order(al[[1]]),]
 ag=ag[order(ag[[1]]),]
 #al[order(ag[[2]], decreasing=T),]
 al=al[order(ag[[2]], decreasing=T),]
 ag=ag[order(ag[[2]], decreasing=T),]
 alxb=dfx[which(dfx$CC.Line %in% al[[1]][which(al[[2]]!=al[[3]])]), c(phen, "Sex")]
 alxc=dfx[which(dfx$CC.Line %in% al[[1]][which(al[[2]]==al[[3]] & al[[2]]=="C")]), c(phen, "Sex")]
 alxa=dfx[which(dfx$CC.Line %in% al[[1]][which(al[[2]]==al[[3]] & al[[2]]=="T")]), c(phen, "Sex")]

 alxa$al="A"
 alxb$al="B"
 alxc$al="C"

 dfz=rbind(alxa,alxb,alxc)
 dfz$at=NA
 dfz$at[which(dfz$Sex=="M" & dfz$al=="A")]=1
 dfz$at[which(dfz$Sex=="F" & dfz$al=="A")]=2
 dfz$at[which(dfz$Sex=="M" & dfz$al=="B")]=3
 dfz$at[which(dfz$Sex=="F" & dfz$al=="B")]=4
 dfz$at[which(dfz$Sex=="M" & dfz$al=="C")]=5
 dfz$at[which(dfz$Sex=="F" & dfz$al=="C")]=6

 write.table(dfz,file=paste("al",phen,"Sx.txt", sep=""),sep="\t",quote=FALSE,row=FALSE)
 return(dfz) 

}


 phen="BMC"
 dfz=alPhen(phen=phen, al=bvs) 
 boxplot(dfz[[1]]~dfz$at, col=c("#006699","#9933CC"), axes=F, outline=F, notch=F, main=phen)
 axis(side=1, at=c(0,1.5,3.5,5.5,7), labels=c("","TT","TC","CC",""))
 stripchart(dfz[[1]]~dfz$at, col=c("blue", "pink"), vertical=T, add=T, method="jitter", jitter=0.2, pch=20, main=phen)


 ph=c("BV.TV", "Tb.NPl", "Tb.Th", "Conn.D", "SMI", "Sp", "Cort.Th", "vBMD", "BMC")
 pdf("~/alDist@chr11.116.pdf") 
 par(mfrow=c(3,3))
 for (phen in ph) { dfz=alPhen(phen=phen, al=bvs);  boxplot(dfz[[1]]~dfz$at, col=c("#006699","#9933CC"), axes=F, outline=F, notch=F, main=phen);  axis(side=1, at=c(0,1.5,3.5,5.5,7), labels=c("","TT","TC","CC",""));  stripchart(dfz[[1]]~dfz$at, col=c("blue", "pink"), vertical=T, add=T, method="jitter", jitter=0.2, pch=20, main=phen)}
 dev.off()
