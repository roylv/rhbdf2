#!/usr/local/bin/Rscript
setwd("/scratch200/roeilevi/rhbdf/")

dt=read.table("http://mus.well.ox.ac.uk/preCC/MEGA_MUGA/Mar2015.MEGA+MDA+MUGA/chr11.MEGA+MDA+MUGA.data", header=F)

n=c(21397:21398,21431:21432)
n2=c(21233:21234,21269:21270 ,n)
n=n2


file="/scratch200/roeilevi/rhbdf/FemOctNoOut.txt"; phen="BV.TV"; thresh=0.17; n=21398; al1="T"; al2="C"; tr=0.0; mid=1



findMtch=function(file, phen, thresh, n, al1, al2, tr, mid=NULL){
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

                        	  if(!is.null(mid)) { 
					if(ag[[2]][i]<(thresh-tr) & dt[nl,-(1:6)][n]!=dt[nl,-(1:6)][n-1]) sd3=c(sd3,(ag[[1]][i]))
					}
			  print(i)
                         }
                 }
       }

        return(list("sd1"=sd1,"sd2"=sd2,"sd3"=sd3))


}
a2=list()
a2=findMtch(file="FemOctNoOut.txt",phen=phen, thresh=thresh, n=n, al1=al1, al2=al2, tr=tr, mid=mid)

save(a2, file="~/findMtch1.RData")

