phens=colnames(df1)[c(10, 12:14, 16,18,24:26,28,29,30,32:37,39:44)]

cors=function(phens, phen2) {
  for(phen in phens){
	ag1=aggregate(df1[[phen]], by=list(df1$CC.Line), FUN=mean, na.rm=T)
	ag2=aggregate(df2[[phen2]], by=list(df2$CC.Line), FUN=mean, na.rm=T)
	 ag2m=ag2[which(ag2[[1]] %in% ag1[[1]]),]
	 ag1m=ag1[which(ag1[[1]] %in% ag2m[[1]]),]


	n1= which(ag1m[[1]] %in% ag1m[[1]][which(is.na(ag1m[[2]]))])
	n2= which(ag2m[[1]] %in% ag1m[[1]][which(is.na(ag1m[[2]]))])
	n3= which(ag1m[[1]] %in% ag2m[[1]][which(is.na(ag2m[[2]]))])
        n4= which(ag2m[[1]] %in% ag2m[[1]][which(is.na(ag2m[[2]]))])

	cat("; n1 = ",n1,"; n2 = ",n2,"; n3 = ", n3,"; n4 = ", n4, "\n")
	
	if(length(n1) >0 & length(n2)>0){
 		ag1m=ag1m[-n1,]
	        ag2m=ag2m[-n2,]
	}	

        if(length(n3) >0 & length(n4)>0){
		if(length(n1)> 0){
			if (sum(n3 %in%  n1)==0) {
        	          ag1m=ag1m[-n3,]
               		  ag2m=ag2m[-n4,]
		      }
      		}
   		else {
		     ag1m=ag1m[-n3,]
                     ag2m=ag2m[-n4,]


		}

	  }





	cr=cor(ag1m[[2]],ag2m[[2]])
	fit=lm(ag2m[[2]]~ag1m[[2]])
	r2=summary(fit)$r.squared
	print(cr)

	pv=anova(fit)$'Pr(>F)'[1]
#pv=lmp(fit)
	mp = mean(scale(ag1m[[2]]))
	hl = mean(scale(ag2m[[2]]))
	lm1=lm(scale(ag2m[[2]])~scale(ag1m[[2]]))[[1]][[1]]
	lm2=lm(scale(ag2m[[2]])~scale(ag1m[[2]]))[[1]][[2]]

#	hl=lm2*mp+lm1


	plot(scale(ag1m[[2]]),scale(ag2m[[2]]), pch=19, lwd=2, axes=F, xlab=phen, ylab=phen2)
	axis(1, at=c(-4:4), labels=c(-4:4))
	axis(2, at=c(-4:4), labels=c(-4:4))
	abline(lm(scale(ag2m[[2]])~scale(ag1m[[2]])), col="blue", lwd=3)
	abline(h=hl, col="black")
	abline(v=mp, col="black")
	legend("topleft",legend=c(paste("r = ",round(cr,3)),
                              paste("r2 = ",round(r2,3)),paste("P = ", round(pv,3))), cex=0.6)

	}
  }

