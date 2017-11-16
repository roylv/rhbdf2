setwd("/scratch200/roeilevi/rhbdf")
source("../merge.19102015.R")
load("mtx2n.RData")
condensed.db = load.condensed.database("../Mar2015.MEGA+MDA+MUGA/CONDENSED/")


for (i in 1:nrow(mtx2)){
	load(paste(mtx2[i,2],mtx2[i,3],"N.RData",sep=""))
	merge.genome(condensed.db=condensed.db, phenotype.name=mtx2[i,2], file=mtx2[i,1], covariate=c("Sex", "WeekOld", "Batch"), regions=T, use.weights=T, mc.cores=20, sb=mtx2[i,3])

}

