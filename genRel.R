map=read.delim("http://mus.well.ox.ac.uk/preCC/MEGA_MUGA/Mar2015.MEGA+MDA+MUGA/chr11.MEGA+MDA+MUGA.map", header=T)
dat=read.table("http://mus.well.ox.ac.uk/preCC/MEGA_MUGA/Mar2015.MEGA+MDA+MUGA/chr11.MEGA+MDA+MUGA.data", header=F)
w=536



file="~/BVTV36675.txt"
d.f=read.delim(file)
sc=read.delim("BV.TV.scan.txt")
sc$marker[which.max(sc$logP)]


#line.df=data.frame(matrix(NA,ncol=1+length(n2),nrow= length(unique(d.f$CC.Line))))

i=1
n=c(21397:21398,21431:21432)
n2=c(21233:21234,21269:21270 ,n)

line.df=data.frame(matrix(NA,ncol=1+length(n2),nrow= length(unique(d.f$CC.Line))))
colnames(line.df)[1:(1+length(n2))]=c("Line",1:length(n2))
i=1
for (line in unique(d.f$CC.Line)) {
	ln=which(dat[,2]==line)	
	line.df$Line[i]=line
	for (j in 2:(length(n2)+1)) line.df[i,j]=as.character(dat[ln,-(1:6)][n2][1,j-1])
        #line.df$a2[i]=as.character(dat[ln,-(1:6)][n][1,2])
        #line.df$b1[i]=as.character(dat[ln,-(1:6)][n][1,3])
        #line.df$b2[i]=as.character(dat[ln,-(1:6)][n][1,4])

	print(line.df[i,])
	i=i+1
}
#aggregate( d.f$BV.TV, by=list(CC.Line=d.f$CC.Line), FUN="mean")
d.f2=d.f[-which(is.na(d.f$BV.TV)),]

ag=aggregate( d.f2$BV.TV, by=list(CC.Line=d.f2$CC.Line), FUN="mean")
#d.f2=d.f[-which(is.na(d.f$BV.TV)),]	
line.df=line.df[order(ag[[2]], decreasing=T),]
ag=ag[order(ag[[2]], decreasing=T),]


ldf2=ldf[-c(1,13,26),]


for (i in 1:31) {
	for(j in 2:5) {
	if (ldf2[i,j]=="C") ldf2[i,j]=1 
	else ldf2[i,j]=2
	
  }
}


lst2=matrix(nrow=nrow(line.df)); for (i in 1:nrow(line.df)) lst[i]=sum(line.df[1,2:9]==line.df[i,2:9])/8




Find matches based on genome and phenoytpe, at a selected locus

findMtch=function(file, phen, thresh, n, al1, al2, tr){
	sd1=NULL
	sd2=NULL
	sd3=NULL
	dfm=read.delim(file)
	dfm=dfm[-which(is.na(dfm[[phen]])),]
	ag=aggregate( dfm[[phen]], by=list(CC.Line=dfm$CC.Line), FUN="mean")
	#print("1")
	for (i in 1:nrow(ag)) {
 		 nl=which(dt[,2]==as.character(ag[[1]][i])) 
		 if (sum(nl)>0) {
			if (!is.na(dt[nl,-(1:6)][n])) {
			  if(ag[[2]][i]>thresh & sum(dt[nl,-(1:6)][(n-1):n]==al1)==2) sd1=c(sd1,ag[[1]][i])
 			  else if (ag[[2]][i]<(thresh-tr) & sum(dt[nl,-(1:6)][(n-1):n]==al2)==2) sd2=c(sd2,ag[[1]][i]) 	
			  else if (ag[[2]][i]<(thresh) & ag[[2]][i]>(thresh-tr) & dt[nl,-(1:6)][n-1]!=dt[nl,-(1:6)][n]) sd3=c(sd3,ag[[1]][i])
			  print(i)
			 }
			}
		  }
	return(list("sd1"=sd1,"sd2"=sd2,"sd3"=sd3))


}
a=list()
a=findMtch(file="FemOctNoOut.txt",phen=phen, thresh=thresh, n=n, al1=al1, al2=al2, tr=tr)

save(a, file="findMtch1.RData")


scan.phenotypes.lm2(smn=5, phn="vBMD",chrMax="chr4", rand=1, condensed.db=condensed.db, data.file="~/BVTV36675.txt", phen="vBMD", permute=150, covariates = c("Sex", "WeekOld","Batch" ), aggregate.data="residuals", use.weights=T, normalise=FALSE, zero.na=FALSE, cousins=TRUE, scan.file="scans_vBMD.RData", mc.cores=40,n=1,mast="FemOctNoOut.txt")
