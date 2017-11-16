sng=read.csv("http://www.sanger.ac.uk/sanger/Mouse_SnpViewer_Export/SNPs/rel-1505?sv=31;sn=2523136;st=68719476735;loc=11:116000000-117000000;format=csv")

mrg=read.delim("BV.TV.chr11:102600000-121700000.merge.txt", sep="")   

#alt way to read sng:
bpr=as.integer(range(mrg$bp[which(mrg$logP.interval>6)]))

loc=paste(mrg$chr[1],":",bpr[1],"-",bpr[2],sep="")
add=paste("http://www.sanger.ac.uk/sanger/Mouse_SnpViewer_Export/SNPs/rel-1505?sv=31;sn=2523136;st=68719476735;loc=",loc,";format=csv", sep="")
sng=read.csv(add)


bpr=range(mrg$bp[which(mrg$logP.interval>6)])
sngx=sng[,c(1:3)]
sngx2=sngx[which(sngx$Position>bpr[1] & sngx$Position<bpr[2]),]

lg=unique(sngx2$Gene)


range(sngx2$Position[which(sngx2$Gene==lg[5])])

lst1=NULL
lst2=NULL
#as.integer(bpr=range(mrg$bp[which(mrg$logP.interval>6)]))
sngx=sng[,c(1:3)]
sngx2=sngx[which(sngx$Position>bpr[1] & sngx$Position<bpr[2]),]
lg=unique(sngx2$Gene)
for (gene in lg) {
   rng=range(sngx2$Position[which(sngx2$Gene==gene)])
   mrgBp=which(mrg$bp>rng[1] & mrg$bp<rng[2])
   prop=sum(mrg[mrgBp,"logP.merge"]>(mrg[mrgBp,"logP.interval"]))/length(mrg[mrgBp,"logP.merge"])
   #print(prop*100) 
   lst1=c(lst1,prop)
   lst2=c(lst2,gene)	
}

gnmtx=matrix(nrow=length(lst1),ncol=2); gnmtx[,1]=lst1; gnmtx[,2]=lst2; gndfbm2=data.frame(gnmtx)

