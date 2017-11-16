#bv=read.delim("scan_1192.txt")
#con=read.delim("scan_ConnD14_408.txt")
#tnh=read.delim("scan_TbTh601.txt")

#wb=1452
#wc=8602
#wn=5123
#wh=4535
#wbl=548
#wnl=154
#whl=185
#wcl=286

setwd("/scratch200/roeilevi/rhbdf")

load("mtx2n.RData")

phen.data=read.delim("~/BVTV36675.txt")
#phens=c("BV.TV", "Tb.NPl", "Conn.D", "Tb.Th", "SMI", "Sp", "vBMD", "Cort.Th")
#phen.data=bv
#phen="BVTV"

source("happyEx19.preCC.R")
condensed.db = load.condensed.database("../Mar2015.MEGA+MDA+MUGA/CONDENSED/")

regHer=function(db, mtx){
#	agg.phen.data = aggregate( phen.data[[phen]],
#                           by=list(CC.Line=phen.data$CC.Line),
#                           FUN="mean" )
	
	hr2=NULL; pv=NULL; phn=NULL; reg=NULL
	for (i in 1:nrow(mtx2)){
		phen.data=read.delim(mtx[i,1])
                agg.phen.data = aggregate( phen.data[[mtx[i,2]]],
          		by=list(CC.Line=phen.data$CC.Line),
                        		   FUN="mean" )

		w=mtx[i,3]
		chr=mtx[i,4]
		int = intersect( db$subjects, unique(agg.phen.data$CC.Line))
		phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]
		db.use = match( phen.data$CC.Line, db$subjects)

		fit.null=lm(x~1,data=phen.data, y=TRUE)
		for (j in 1:length(db$additive$chr[[chr]])) {
			e=db$additive$chr[[chr]][[as.numeric(j)]]$mat
			d=e[db.use,]	
			fit.genetic=lm(x~d, data=phen.data)
			a=anova(fit.null, fit.genetic)
			hr2=c(hr2,(a[1,2]-a[2,2])/a[1,2]); pv=c(pv,a[2,6]); phn=c(phn, mtx[i,2])
			cat("Hr2=",(a[1,2]-a[2,2])/a[1,2], "Pv=",a[2,6], "\n", sep="\t")
		}
		
	}
	lst=list("phen"=phn,"hr2"=hr2,"pv"=pv); save(lst,file="regh2resLoc.RData")
}


regHer(db=condensed.db, mtx=mtx2)


