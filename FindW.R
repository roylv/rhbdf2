setwd("/scratch200/roeilevi/rhbdf")
load("mtx2n.RData")
source("happyEx19.preCC.R")

condensed.db = load.condensed.database("../Mar2015.MEGA+MDA+MUGA/CONDENSED/")

w=NULL
chr=NULL
lp=NULL
for (i in 1:nrow(mtx2)){
	s=scan.phenotypes.lm(condensed.db, data.file=mtx2[i,1], phen=mtx2[i,2], permute=200, histogram.pdf=NULL, covariates = c("Sex", "WeekOld","Batch"), mc.cores=40, aggregate.data="residuals",use.weights=T, scan.file=paste(mtx2[i,2],mtx2[i,3], "N.RData",sep=""))
	scan=s[[mtx2[i,2]]]$result
	chrx=scan[which.max(scan$logP),"chromosome"]	
	chr=c(chr,chrx)
	wx=which.max(scan[which(scan$chromosome==chrx),"logP"])
	w=c(w,wx)
	lp=c(lp, max(scan$logP))
	}
save(w,file="loci.RData")
save(chr, file="chrs.RData")
save(lp, file="logPs.RData")
