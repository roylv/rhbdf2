library(happy.hbrem)
#library(multicore)
library(parallel)

library(survival)
library(hash)

# HMM FUNCTIONS

make.condensed.genome <- function( chrs, prefix, gdir= "./", dir="./CONDENSED/", thin=20, extension=".data", gen=7, haploid=FALSE, mc.cores=10) {

    if ( mc.cores == 1 )
        lapply( chrs, make.condensed.chr.database, prefix=prefix, gdir=gdir, dir=dir, thin=thin, extension=extension, gen=gen, haploid=haploid)
    else
        mclapply( chrs, make.condensed.chr.database, prefix=prefix, gdir=gdir, dir=dir, thin=thin, extension=extension, gen=gen, haploid=haploid, mc.cores=mc.cores)
  if ( haploid ) 
    models = c( "additive" )
  else 
    models = c( "additive", "full" )

  db = list( chrs=chrs, prefix=prefix, dir=dir, thin=thin, extension=extension, gen=gen, haploid=haploid, models=models, chrs=chrs )
  save( db, file=paste( dir, "/db.RData", sep=""))
}


make.condensed.chr.database <- function( chr, prefix, gdir="./", dir="./CONDENSED/", thin=20, extension=".data", gen=7, haploid=FALSE) {

  dir.create(dir, recursive=TRUE, show=FALSE) 
  
  pedfile = paste(gdir, chr, prefix, extension, sep="")
  allelesfile = paste( gdir, chr, prefix, ".alleles", sep="")
  mapfile = paste( gdir, chr, prefix, ".map", sep="")
  cat( pedfile, file.exists(pedfile), allelesfile, file.exists(allelesfile), mapfile, file.exists(mapfile), "\n")
  h = happy( pedfile, allelesfile, gen=gen, file="ped", mapfile=mapfile, haploid=haploid)

  if ( haploid ) 
    models = c( "additive" )
  else 
    models = c( "additive", "full" )

  subjects=h$subjects
  strains=h$strains
  save(subjects, file=paste(dir, "/subjects.RData", sep=""))	
  save(strains, file=paste(dir, "/strains.RData", sep=""))	
  condensed.params = c(thin=thin, pedfile=pedfile, allelesfile=allelesfile, mapfile=mapfile, gen=gen, haploid=haploid, models=models)
  save(condensed.params, file=paste(dir, "/", chr, ".param.RData", sep=""))		
  subdir = paste(dir, "/chr/", sep="")
  dir.create(subdir, recursive=TRUE, show=FALSE) 

  for( model in models ) {
    scandir = paste(subdir, model, "/", sep="")	
    dir.create(scandir, recursive=TRUE, show=FALSE) 

    condensed = make.condensed.matrices( h, thin, model=model )
    Rfile = paste(scandir, chr, ".RData", sep="")  
    save(condensed, file=Rfile)
    x = sapply( condensed, function( x ) { return (c( x$from, x$to, x$from.marker, x$to.marker, x$len ))})
    condensed.summary = data.frame(t(x))
    names(condensed.summary) = c( "from", "to", "from.marker", "to.marker", "width")

    condensed.summary$from.bp=h$bp[match(condensed.summary$from.marker,h$marker)]
    condensed.summary$to.bp=h$bp[match(condensed.summary$to.marker,h$marker)]
    condensed.summary$chromosome=rep(chr,nrow(condensed.summary))
    save(condensed.summary, file=paste(scandir, chr, ".condensed.summary.RData", sep=""))
    genome = data.frame(marker=h$markers,chromosome=h$chromosome,bp=h$bp,cm=h$map)
    save(genome, file=paste(scandir,  chr, ".summary.RData", sep=""))
  }
  
}	

                                        # loads a pre-created database of condensed HAPPY design matrices into memory
                                        # returns an object of class "condensed.happy" with the following fields
                                        # subjects
                                        # strains
                                        # chr$additive$$chr$chr1 ....chr$additive$chr19 ....
                                        # additive$condensed.summary
                                        # chr$full$chr$chr1 ....chr$full$chr19 ....



load.condensed.database <- function( dir="./CONDENSED.20/" ) {

  load( paste(dir, "/dbEx19.RData",sep=""))

  load( paste(dir, "/strains.RData", sep=""))
  db$strains = strains;
  
  load( paste(dir, "/subjects.RData", sep=""))
  db$subjects = subjects;

  cat( "loading condensed db", dir, "\n")

  models = db$models 
  for( model in models ) {
    n = 0;
    db[[model]] = list()
    db[[model]]$chr = list()
    db[[model]]$matrices = list()
    h = hash()
    files = paste(dir, "/chr/", model, "/", db$chrs, ".RData", sep="")
    for ( i in 1:length(files)) {
      chr.data = files[i]
      load( chr.data )
      db[[model]]$chr[[db$chrs[i]]]=condensed
      db[[model]]$matrices = c( db[[model]]$matrices, condensed) 
      n = n + length(db[[model]]$chr[[db$chrs[i]]])
      lapply( condensed, function( item ) { h[[item$from.marker]] <<- item$mat } )
    }
    db[[model]]$hash = h
    cat( "model", model, "read ", n, " matrices\n");
    cat( "loading genome summary\n")
    db[[model]]$summary = NULL
    files = paste(dir, "/chr/", model, "/", db$chrs, ".summary.RData", sep="")
    for (chr.summary in files ) {
      load( chr.summary )
      db[[model]]$summary = rbind( db[[model]]$summary, genome )
    }

    cat ( "loading condensed summary\n")
    db[[model]]$condensed.summary = NULL
    files = paste(dir, "/chr/", model, "/", db$chrs, ".condensed.summary.RData", sep="")
    for (condensed.chr.summary in files ) {
      load( condensed.chr.summary )
      db[[model]]$condensed.summary = rbind( db[[model]]$condensed.summary, condensed.summary )
    }

# create some data objects so that the resulting object looks like a genome.cache...

    idx = match( db[[model]]$condensed.summary$from.marker, db[[model]]$summary$marker)
    db[[model]]$subjects = subjects
    db[[model]]$subjects = strains
    db[[model]]$genome = db[[model]]$summary[idx,]
    db[[model]]$genome$map = db[[model]]$genome$cm
    db[[model]]$genome$chr = gsub( "chr", "", db[[model]]$genome$chromosome )
    db[[model]]$markers = db[[model]]$genome$marker
    db[[model]]$chromosome = db[[model]]$genome$chromosome
    db[[model]]$map = db[[model]]$genome$map
  }
  load( paste(dir, "/strains.RData", sep=""))
  db$strains = strains;
  
  load( paste(dir, "/subjects.RData", sep=""))
  db$subjects = subjects;
  
  db$markers = db$additive$markers
  class(db) = "condensed.happy";
  return(db)
}

condensed.hdesign <- function( db, marker, model="additive" ) {
	if ( class(db) == "condensed.happy" ) 
		return(db[[model]]$hash[[marker]])
	else 
		return( hdesign( db, marker, model ))
}

# LINEAR MODELS

scan.phenotypes.lm <- function( condensed.db = NULL, data.file="NA", phen, permute=0, histogram.pdf=NULL, covariates = c("Sex"), aggregate.data="residuals", use.weights=FALSE, zero.na=FALSE, cousins=TRUE, scan.file="scans.RData", mc.cores=10, SEX=SEX ) {

    if ( is.null(condensed.db)) condensed.db = load.condensed.database("/Net/mus/data/www/preCC/MEGA_MUGA/Feb2014.MEGA+MDA+MUGA/CONDENSED/")
    
#   if (!is.null(SEX)){

	#if(SEX=="M"){
	#phen.data = read.delim(data.file)
	#phen.data = subset(phen.data, Sex=="M")
	#phen.data$Sex = factor(phen.data$Sex)
	#print(phen.data$Sex)
	#}
	#else if (SEX=="F" ){
 	#phen.data = read.delim(data.file)
        #phen.data = subset(phen.data, Sex=="F")
        #phen.data$Sex = factor(phen.data$Sex)
	#print(phen.data$Sex)
#		}
#	}
	
#    else {
	
	phen.data = read.delim(data.file)
	#phen.data<-f
#   	}	
    err=0
    include.sex = FALSE

    cols = c( phen, covariates )
    for ( p in cols ) {
        if ( is.null( phen.data[[p]] ) ) {
            warning( "file ", data.file, " has no column ", p, "\n")
            err = err+1
            if ( p == "Sex" ) include.sex=TRUE
        }
    }
    
    if ( err > 0)  return(NULL)
    if ( "aggregate.traits" == "sex" & ! include.sex ) {
        warning(" error , no sex\n")
        return(NULL)
    }

    if ( aggregate.data == "residuals") {
	null.model = "~ 1"
    }
    else if ( aggregate.data == "sex") {
	null.model = "~ Sex"
    }
    
    if ( cousins==FALSE ) {
	cous = cousin.lines( phen.data$CC.Line, rm=TRUE )
	dm1 = nrow(phen.data)
        browser()
        phen.data = phen.data[match(cous$id,phen.data$CC.Line),]
	dm2 = nrow(phen.data)
	cat( "removing cousins from", dm1, "to", dm2, "\n")
    }
    
    scans = list()
 
    phens = c()
    agg.data = list()
    for ( p in phen ) {
        cat( p, "\n")
        
	if ( is.numeric(phen.data[[p]]) ) {
            phens = c(phens, p)
            if ( zero.na==TRUE ) {
                cat("setting 0 to NA\n")
                n0 = sum(phen.data[[p]] == 0.0)
                phen.data[[p]] = ifelse( phen.data[[p]] == 0.0 , NA, phen.data[[p]] )
                na = sum(is.na(phen.data[[p]]) )
                cat("setting 0 to NA", n0, na, "\n")
            }
          
	    if ( aggregate.data == "residuals"  ) {
                agg.phen.data = aggregate.phenotype(data.file, p, covariates, SEX=SEX)
                agg.data[[p]] = agg.phen.data

            }
            else if ( aggregate.data == "sex" ) {
                complete.sex = lm( phen.data[[p]] ~ Sex, data=phen.data)
                complete.sex.line = lm( phen.data[[p]] ~ Sex + CC.Line, data=phen.data )
                cat( "complete ANOVA for ", p,"\n")
                print(anova(complete.null,complete.sex,complete.sex.line))
          
                agg.phen.data = aggregate( phen.data[[p]], by=list(CC.Line=phen.data$CC.Line, Sex=phen.data$Sex), FUN="mean")#FUN="median")
                agg.phen.data[[p]] = agg.phen.data$x
                agg.data[[p]] = agg.phen.data
                
            }
         
	  else { # leave input phenotype data unchanged
                agg.phen.data = phen.data
                agg.data = agg.phen.data
            }
         
	   # if (!use.weights) {	
           # complete.null = lm( agg.phen.data[[p]] ~ 1 )
	   # }
	   # else {#(use.weights==TRUE){
           # complete.null = lm( agg.phen.data[[p]] ~ 1 , weights=agg.data[[p]]$weights) #  wgg.phen.datagg.phen.dataeights=agg.data[[p]]$weights )

	#	print(agg.data[[p]]$weights)
           # }
#      scan = condensed.genome.scan.lm( condensed.db, agg.phen.data, paste( p, null.model), permute=permute, phen.name=phen.name)
            scan = multicore.condensed.genome.scan.lm( condensed.db, agg.phen.data, paste( p, null.model), permute=permute, mc.cores=mc.cores, use.weights=use.weights)
            write.table(scan$result, file=paste(p, ".scan.txt", sep=""), quote=FALSE, row=FALSE, sep="\t")
            scans[[p]] = scan
        }
        else {
            warning( "error phenotype ", p, "not numeric\n")
        }
    }
  
    plot.histograms( histogram.pdf, agg.phen.data, phens  ) 
    if ( !is.null(scan.file ) ) {
        cat( "writing scans to ", scan.file, "\n")
        scan.data=list( scans=scans, data.file=data.file, phen=phen, permute=permute,  covariates=covariates, aggregate.data=aggregate.data, agg.data=agg.data )
        save( scan.data, file=scan.file )
    }
    
    return(scans)
}


condensed.genome.scan.lm <- function( condensed.db, phen.data, null.model, permute=0, model="additive", phen.name = NULL ) {

  if ( class(condensed.db) != "condensed.happy" ) {
    warning( "condensed.db not of class condensed.happy\n")
    return(NULL)
  }


  int = intersect( condensed.db$subjects, unique(phen.data$CC.Line))
  cat( length(int), "subjects analysed with phenotypes\n")


#  phen.data = phen.data[phen.data$CC.Line %in% int,]
  phen.data = phen.data[match( int, phen.data$CC.Line ,nomatch=0),]

  if ( !is.null(phen.name)) 
    hist(phen.data[[phen.name]], xlab=phen.name, main=phen.name)

  permuted.data=NULL
  if ( is.numeric(permute) & permute>0) {
    f = lm( as.formula(null.model), data=phen.data, y=TRUE)
    permuted.data = replicate( permute, sample(residuals(f)))
    permuted.data = cbind( residuals(f), permuted.data)
  }
  
#  db.use = match(int, phen.data$CC.Line, nomatch=0)
  db.use = match( phen.data$CC.Line, condensed.db$subjects)

  fit.null = lm( null.model, data=phen.data );
  print(anova(fit.null))

  genetic.model = paste( as.character(null.model), " + d")
  result = condensed.db[[model]]$condensed.summary


  if ( is.numeric(permute) & permute>0) {
    chr = condensed.db[[model]]$chr
    perm.result = NULL

    for ( i in 1:length(chr)) {
      mats = chr[[i]]
      plp = scan.chromosome.lm( mats, db.use, genetic.model, phen.data, permuted.data, fit.null, condensed.db$subjects)
      perm.result = rbind(perm.result,t(plp))
    }

    result$logP = perm.result[,1]
    
    result$pointwise.permute.pval = apply( perm.result, 1, function(w) {sum(w[1]<=w[2:length(w)])} )/permute
    genomewide.distribution = sort(apply( perm.result[,2:ncol(perm.result)], 2, max ))
    result$genomewide.permute.pval = sapply( result$logP, function(x,genomewide.max) { sum(x <= genomewide.max)}, genomewide.max)/permute
    return( list(result=result, genomewide.distribution=genomwide.distribution) )
  }
  else {
    logP = NULL
    chr = condensed.db[[model]]$chr

    for ( i in 1:length(chr)) {
      mats = chr[[i]]
      lp = sapply( mats, function( e, db.use, genetic.model, phen.data, fit.null, subjects) {

        d = e$mat
        rownames(d) = subjects
        d = d[db.use,]
        fit.genetic = lm( as.formula(genetic.model), data=phen.data) 
        a = anova(fit.null, fit.genetic)
        return( -log10(a[2,6]))
      }
        , db.use, genetic.model, phen.data, fit.null, condensed.db$subjects)
      logP = c( logP, lp )		
    }
    result$logP = logP
    return(list(result=result))
  }  

}

scan.chromosome.lm <- function( mats, db.use, genetic.model, phen.data, permuted.data, fit.null, subjects, use.weights=use.weights ) {
  return( sapply( mats, test.locus.lm, db.use, genetic.model, phen.data, permuted.data, fit.null, subjects, use.weights=use.weights ))
}




test.locus.lm <- function( e, db.use, genetic.model, phen.data, permuted.data, fit.null, subjects, use.weights=use.weights) {

#	browser()
  d = e$mat
  rownames(d) = subjects
  d = d[db.use,]
	if (!use.weights) {
  fit.genetic = lm( as.formula(genetic.model), data=phen.data, qr=TRUE, y=TRUE)
  }
else {
	fit.genetic = lm( as.formula(genetic.model), data=phen.data, qr=TRUE,y=TRUE, weights=weights)
}
  a = anova(fit.null, fit.genetic)
  y = fit.genetic$y
  y = y-mean(y)
  tss = sum(y*y)
  n = length(y)-fit.null$rank
  df = fit.genetic$rank - fit.null$rank
  rss = apply( permuted.data, 2, function( y, qr ) { r = qr.resid(qr,y); return(sum(r*r)) } , fit.genetic$qr )
  fss = tss-rss; 
  f.stat=(fss/df)/(rss/(n-df)); 
  
  p = -pf( f.stat, df, n-df, lower=FALSE, log=TRUE)/log(10)
 # print("DEBUGGINH\n")
#	print(a)
#print("genetic")
#print(anova(fit.genetic))
#print("null")
#print(anova(fit.null))
#	print(tss)
#	print(fss)
#print(f.stat)
#print(f.stat)
#save.data = list( phen.data=phen.data, d=d)
#save(save.data,file="save.RData")

	#print ("END") 
 return( c(-log10(a[2,6]), p))
 return( p)
}



multicore.condensed.genome.scan.lm <- function( condensed.db, phen.data, null.model, permute=0, model="additive", mc.cores=8, subject.name="CC.Line", use.weights=use.weights) {

  if ( class(condensed.db) != "condensed.happy" ) {
    warning( "condensed.db not of class condensed.happy\n")
    return(NULL)
  }
  

  int = intersect( condensed.db$subjects, unique(phen.data[[subject.name]]))
  cat( length(int), "subjects analysed with phenotypes\n")

  #print(phen.data)
  phen.data = phen.data[match(int, phen.data[[subject.name]],nomatch=0),]
  #phen.data=as.table(phen.data, sep=",")
  permuted.data=NULL
  if ( is.numeric(permute) & permute>0) {
    f = lm( as.formula(null.model), data=phen.data, y=TRUE)
    permuted.data = replicate( permute, sample(residuals(f)))
    permuted.data = cbind( residuals(f), permuted.data)
  }
  #phen.data$weights<-as.numeric(phen.data$weights) 
  #w<-phen.data$weights
  db.use = match(int, condensed.db$subjects,nomatch=0)
  #w<-phen.data$weights
  #phen.data[,1] = as.numeric(phen.data[,1])
  if( !use.weights) {
 
	fit.null = lm( null.model, data=phen.data )
	}
  else {#if (!is.null(phen.data$weights)) {
	
	#d1=as.data.frame(phen.data)
	#print(phen.data$weights)
	#w=as.numeric(w)
	#print(w)
	#w=as.numeric(as.character(sub("," , ".", phen.data$weights)))
	#print(w)
        #print(str(phen.data))
        #print(null.model)
	fit.null = lm( null.model, data=phen.data, weights=weights)
	#print("XXXXXXXXXXX")
	
	}   
  
  print(anova(fit.null))

  genetic.model = paste( paste(as.character(null.model)), " + d", collapse=" ")
  result = condensed.db[[model]]$condensed.summary

  if ( is.numeric(permute) & permute>0) {
    chr = condensed.db[[model]]$chr
    perm.result = NULL
    if ( mc.cores > 1 ){ 
	    plps = mclapply(chr, scan.chromosome.lm, db.use, genetic.model, phen.data, permuted.data, fit.null, condensed.db$subjects, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores=mc.cores, use.weights=use.weights )
	}
    
    else { 
	plps = lapply(chr, scan.chromosome.lm, db.use, genetic.model, phen.data, permuted.data, fit.null, condensed.db$subjects )
	}
    for ( i in 1:length(plps)) {
      perm.result = rbind(perm.result,t(plps[[i]]))
    }

    result$logP = perm.result[,1]
    
    result$pointwise.permute.pval = apply( perm.result, 1, function(w) {sum(w[1]<=w[2:length(w)])} )/permute
    genomewide.distribution = sort(apply( perm.result[,2:ncol(perm.result)], 2, max ))
    result$genomewide.permute.pval = sapply( result$logP, function(x,genomewide.distribution) { sum(x <= genomewide.distribution)}, genomewide.distribution)/permute
    return( list(result=result, genomewide.distribution=genomewide.distribution) )
  }
  else {
    logP = NULL
    chr = condensed.db[[model]]$chr

    for ( i in 1:length(chr)) {
      mats = chr[[i]]
      lp = sapply( mats, function( e, db.use, genetic.model, phen.data, fit.null, subjects) {

        d = e$mat
        rownames(d) = subjects
        d = d[db.use,]
	if ( is.null(phen.data$weights)) {
        	fit.genetic = lm( as.formula(genetic.model), data=phen.data) 
	}       
	else {
		fit.genetic = lm( as.formula(genetic.model), data=phen.data, weights=weights)
}

a = anova(fit.null, fit.genetic)
        return( -log10(a[2,6]))
      }
        , db.use, genetic.model, phen.data, fit.null, condensed.db$subjects)
      logP = c( logP, lp )		
    }
    result$logP = logP
    return(list(result=result))
  }  

}


# GLM



make.condensed.matrices <-function( h, thin=100, model="additive" ) {
  nm = length(h$markers)-thin
  condensed = lapply( seq(1,nm,thin), function ( m, h, thin ) {
    mat = hdesign( h, m, model=model )
    m2=m
    if ( thin > 1 ) {
      m1 = m+1
      m2 = m+thin-1
      if ( m2 >= length(h$markers)) m2 = length(h$markers)-1
      denom = m2-m+1

      for( i in m1:m2) {
        mat = mat + hdesign( h, i, model=model )
      }
      mat = mat/thin
    }      
                                  #		cat( m, m2, h$markers[m],h$markers[m2],thin, "\n")
    return( list( mat=mat, from=m, to=m2, from.marker=h$markers[m], to.marker=h$markers[m2], len=thin ))
  }, h, thin )
  return(condensed)
}


#PLOTTING

plot.histograms <- function( histogram.pdf, phen.data, phen ) {
  if ( ! is.null(histogram.pdf)) {
    pdf(histogram.pdf)
    par(mfrow=c(3,3))
    for ( p in phen ) {
      if ( ! is.null(phen.data[[p]]))
        hist(as.numeric(as.character(phen.data[[p]])), xlab=p, main=p)
	#hist(as.numeric(as.character(df[[p]])), xlab=p, main=p)
	#hist(as.numeric(phen.data$p$result$logP), xlab=p, main=p)
	#hist(phen.data$p$result, xlab=p, main=p)
    }
    dev.off()
  }
}

plot.qq <- function( qqplot.pdf, phen.data, phen ) {
  if ( ! is.null(qqplot.pdf)) {
    pdf(qqplot.pdf)
    par(mfrow=c(3,3))
    for ( p in phen ) {
      if ( ! is.null(phen.data[[p]]))
        qqnorm(as.numeric(as.character(phen.data[[p]])), xlab=p, main=p)
	qqline(as.numeric(as.character(phen.data[[p]])), xlab=p, main=p)
        #hist(phen.data$p$result, xlab=p, main=p)
    }
    dev.off()
  }
}

plot.box <- function( boxplot.pdf, phen.data, phen ) {
  if ( ! is.null(boxplot.pdf)) {
    pdf(boxplot.pdf)
    par(mfrow=c(3,3))
    for ( p in phen ) {
      if ( ! is.null(phen.data[[p]]))
        boxplot(as.numeric(as.character(phen.data[[p]])), xlab=p, main=p)
        
        #hist(phen.data$p$result, xlab=p, main=p)
    }
    dev.off()
  }
}



#plot.scans <- function ( db, scans, pdffile=NULL,q = c( 0.5, 0.9, 0.95 ) ) {
plot.scans <- function ( scans, pdffile=NULL,q = c( 0.5, 0.9, 0.95 ) ) {

  if ( ! is.null(pdffile)) pdf( pdffile ) 
  par(mfrow=c(4,1),cex=0.4)	
  #DF = data.frame(1:3)
  DF = data.frame(a=character(),b=character(),c=character(),d=character())
  #names(DF)<-c("a","b","c")
  for( p in names(scans)) {
    scan = scans[[p]]$result

    quantiles = scans[[p]]$genomewide.distribution[ floor( length(scans[[p]]$genomewide.distribution)*q) ]
    #print(quantiles)
    #DF = rbind(DF,c([[p]],as.list(quantiles)))
    my = max(quantiles,scan$logP)	
    chr = scan$chromosome[2:nrow(scan)]
    boundary = c(1, which(chr != scan$chromosome[1:nrow(scan)-1]))
    boundary2 = c(boundary[2:length(boundary)],nrow(scan))
    xc = 0.5*(boundary+boundary2)
    if (! is.null(pdffile) ) {
	plot( 1:nrow(scan), scan$logP, t="l", ylim=c(0,my*1.05),main=p, ylab="logP", xaxs="i", xaxt="n", xlab="position")
    	abline( v=boundary, col="red")
    	abline(h=quantiles,col="grey")
    	axis(1, at=xc, labels=chr[boundary])
    }
    if ( !is.null(scans[[p]]$qtls) ) {
	qtls = scans[[p]]$qtls[scans[[p]]$qtls$phenotype==p,]
	x = match( qtls$peak.SNP, scan$from.marker)
	points(x,rep(my,length(x)),pch=20,col="orange")
    }
  }
  if ( ! is.null(pdffile)) dev.off()
  #colnames(DF)<-c("a","b","c","d")
  #return(quantiles)
}


plot.scans2 <- function (s = NULL, phen, q= c( 0.5, 0.9, 0.95 ) ) {
  
  #if ( ! is.null(pdffile)) pdf( pdffile )
  #par(mfrow=c(4,1),cex=0.4)
  if (is.null(s)) load("scans.RData")
  #ph=scan.data$scan
  #scans=paste(ph,"[",phen,"]", sep="")
  scans=scan.data$scan[[phen]]
  #print(names(scans))		
  #rint(scans)
  #DF = data.frame(1:3)
  #DF = data.frame(a=character(),b=character(),c=character(),d=character())
  #names(DF)<-c("a","b","c")
  #for( p in names(scans)) {
  scan = scans[[1]]
  quantiles = scans[[2]][ floor( length(scans[[2]])*q) ]
    #print(quantiles)
    #DF = rbind(DF,c([[p]],as.list(quantiles)))
    #my = max(quantiles,scan$logP)
    #chr = scan$chromosome[2:nrow(scan)]
    #boundary = c(1, which(chr != scan$chromosome[1:nrow(scan)-1]))
    #boundary2 = c(boundary[2:length(boundary)],nrow(scan))
    #xc = 0.5*(boundary+boundary2)
    #if (! is.null(pdffile) ) {
    #    plot( 1:nrow(scan), scan$logP, t="l", ylim=c(0,my*1.05),main=p, ylab="logP", xaxs="i", xaxt="n", xlab="position")
    #    abline( v=boundary, col="red")
    #    abline(h=quantiles,col="grey")
    #    axis(1, at=xc, labels=chr[boundary])
    #}
    #if ( !is.null(scans$qtls) ) {
    #    qtls = scans[[p]]$qtls[scans[[p]]$qtls$phenotype==p,]
    #    x = match( qtls$peak.SNP, scan$from.marker)
    #    points(x,rep(my,length(x)),pch=20,col="orange")
    #}
  #}
  #if ( ! is.null(pdffile)) dev.off()
  #colnames(DF)<-c("a","b","c","d")
  return(quantiles)
}



# DATA IO

write.gscandb.scans <- function( scans, phenotype, sh.file=NULL ) {

        if ( is.null(sh.file))
		sh.file=paste(phenotype,".gscan.upload.sh",sep="")
        handle = file(sh.file,"w")
        cat("cd /data/www/cgi/gscan/\n", file=handle)

	lapply( names(scans), function( p, scans,handle ) {
                cwd = getwd()
		gscandb = data.frame(marker=scans[[p]]$result$from.marker)
		nam = paste( phenotype, ".", p , sep="")	
	       	gscandb$additive.logP = scans[[p]]$result$logP
		gscandb.file = paste(cwd, "/", nam, ".gscandb.txt", sep="")
		cat( "writing ", gscandb.file, "\n")

		cat( "sh gscanDelete.sh population=preCC phenotype=",nam," label=additive.logP\n",sep="", file=handle)
        	cat( "sh gscanLoad.sh pop=preCC build=37 pheno=",nam," labels=additive.logP update gscan=", gscandb.file, "\n", sep="", file=handle)
		write.table(gscandb,file=gscandb.file,sep="\t",quote=FALSE,row=FALSE)
		q = c( 0.5, 0.9, 0.95 ) 
		quantiles = NULL
		if ( !is.null(scans[[p]]$genomewide.distribution) )
		        quantiles = scans[[p]]$genomewide.distribution[ floor( length(scans[[p]]$genomewide.distribution)*q) ]
		else if ( !is.null(scans[[p]]$permuted.logP ) ) 
		        quantiles = scans[[p]]$permuted.logP[ floor( length(scans[[p]]$permuted.logP)*q) ]
		cat( "preCC", "37", nam, "additive.logP", "Q=0.50", quantiles[1], "Q=0.90", quantiles[2], "Q=0.95", quantiles[3], "\n", sep="\t")
        }, scans, handle )

	close(handle)
}

aggregate.phenotype <- function( data.file, phenotype, covariates, FUN="mean", SEX=SEX) {#FUN="media" ) {
    
     
     d = read.delim( data.file )
 #    d = na.omit(d)
 #    if(!is.null(SEX)){
 #    if (SEX=="M"){
	
  #      d = subset(d, Sex=="M")
#	d$Sex = factor(d$Sex)
  #      print(d$Sex)
  #      }
 #   else if (SEX=="F" ){
   #     d=d[d$sex == "F",]
        #d3 = d$Sex
	#print(colnames(d3))
	#print("XXXXXXXXXXXXXXXXXXXXXX")
	#d3 = subset(d3, Sex="F")
	#d$Sex = factor(d$Sex)
	#print(d3)
	#d3 = as.data.frame(d3)
	#print("XXXXXXXX")
	#d = merge(d, d3, by="Sex")
	#colnames(d3) <-c("Sex")
	#d = merge (d, d3, by="Sex")
        #print(d$Sex)
#	}
 #   }
#
    
    d = d[,c(phenotype, covariates, "CC.Line")]
    d = d[complete.cases(d),]
    #print("XXXXXXXXXX")
    #d3 = d[d$Sex=="F",] 
    #print(d3)
    #print("XXXXXXXXX")
    if ( length(covariates) > 0 ) {
        if ( length(covariates) > 1 ) {
            cov = paste( covariates , collapse=" + " )
        }
        else {
            cov = covariates
        }
       
	formula = as.formula( paste( phenotype, "~" , cov))
        #d = as.data.frame(d)
	#print(d$CC.Line)
	#print(is.data.frame(d))
	#print(levels(d))
	#print("XXXXXXXXX")
	#d=subset(d, Sex=="F")
	# d3<- as.data.frame(lapply(d3, function(d3){ordered(d3,
        #      levels = letters[5:1],labels=letters[5:1])}))
	
        #require(plyr)
        #d3<- catcolwise( function(v) ordered(v, levels = letters[5:1]))(d3)
	#print(d3)
	f = lm ( formula, data=d)
	#require(plyr)
	#d3<- catcolwise( function(v) ordered(v, levels = letters[5:1]))(d3)
        #Y <- catcolwise( function(v) ordered(v, levels = letters[5:1]))(X)
	
	#d3<- as.data.frame(lapply(d3, function(d3){ordered(d3,
        #      levels = letters[5:1],labels=letters[5:1])}))
     #   print("XXXXXXXX") 
	cat( "aggregating", phenotype, nrow(d), "observations\n");
        a = aggregate( resid(f), by=list(CC.Line=d$CC.Line), FUN="mean") # change to "mean"
       
        cat(nrow(a), "CC.Line\n");
        print(anova(f))
	#print("XXXXXX")
    }
    else {
        a = aggregate( d[[phenotype]], by=list(CC.Line=d$CC.Line), FUN=FUN)
    }
    # need to count the number of occurences within each line , say in column N
    d2 = subset(d, select=c("CC.Line"))
    #d2 = subset(read.delim(data.file), select=c("CC.Line"))
    d2 = as.data.frame(table(d2))
        
    colnames(d2)<-c("CC.Line", "weights")
    d2$weights=1/d2$weights
    #colnames(d2)[2] <- "weights"    
    #colnames(d2)[1] <- "CC.Line"
       
    
   
    ag = data.frame( CC.Line=a$CC.Line)
    ag[[phenotype]] = a$x
    ag = merge(ag, d2, by="CC.Line")

    #ag$weight = weights

    return(ag)
}

read.phenotypes <- function( data.file, condense=TRUE, condense.phenotype=NULL ) {

  df = read.delim(data.file)
  nam = names(df)
  err=0
  if ( is.null(df$CC.Line) ) {
    warning( "file ", data.file, "has no column CC.Line\n")
    err = err+1
  }
  if ( is.null(df$sex) ) {
    warning( "file ", data.file, "has no column sex\n")
    err = err+1
  }
  if ( err > 0)  return(NULL)
  
  
  if ( condense ) {
    condensed = list()
    CC.lines = unique(as.character(df$CC.Line))

    if ( is.null(condense.phenotype)) {
      for ( col in names(df)) {
        if ( is.numeric(df[[col]]) ) {
	   complete.null = lm( df[[col]] ~ 1 )
   	   complete.sex = lm( df[[col]] ~ Sex, data=df)
   	   complete.sex.line = lm( df[[col]] ~ Sex + CC.Line, data=df )
	   cat( "ANOVA for ", col,"\n")
	   print(anova(complete.null,complete.sex,complete.sex.line))
	

          a = aggregate( df[[col]], by=list(CC.Line=df$CC.Line, sex=df$sex), FUN="mean")#FUN="median")
          condensed[[col]] = a$x[match(CC.lines,a$CC.Line)];
        } else condensed[[col]] = df[match(CC.lines, df$CC.Line),col] 
      }
    }
    else {
      if ( is.numeric(df[[condense.phenotype]]) ) {
        a = aggregate( df[[condense.phenotype]], by=list(CC.Line=df$CC.Line), FUN=function(x) { quantile( x, p=c(0.5), type=1 ) } )
        idx = which( a$CC.Line==a$CC.Line & a[[condense.phenotype]]==a$x)
        
        condensed = a[idx,]
      } else {
        warning( "condensed phenotype ", condense.phenotype, "not numeric\n")
        return(NULL)
      }
    }
    
    return(data.frame(condensed))
  } else {
    return(df)
  }
}




# HBREM

scan.chromosome.hbrem <- function( loci, db.use, phen.data, HaploidInd, Ndip, Nind, Npost, Nbin) {

  return(sapply( loci, test.locus.hbrem, phen.data, db.use, HaploidInd, Ndip, Nind, Npost, Nbin) )
}

test.locus.hbrem <- function(locus, phen.data, db.use, HaploidInd, Ndip, Nind, Npost, Nbin) {

    d = locus$mat[db.use,]

    hb <- hbrem(RX=d, HaploidInd=HaploidInd, Ndip=Ndip, Nind=Nind, Npost=Npost, Nbin=Nbin, Ry=phen.data)

    reg.full.lm <- lm(phen.data ~ d)
    reg.null.lm <- lm(phen.data ~ 1)
    Ftest <- anova(reg.null.lm, reg.full.lm)
    pval <- Ftest$"Pr(>F)"[2]

#    cat( locus$from.marker, -log10(pval), hb[[1]][12],"\n" )
    return( c(  -log10(pval), hb[[1]], hb[[2]], hb[[3]], hb[[4]] ))
  }

multicore.condensed.genome.scan.hbrem <- function( condensed.db, phen.data, mc.cores=8, haploid=FALSE  ) {


  if ( class(condensed.db) != "condensed.happy" ) {
    warning( "condensed.db not of class condensed.happy\n")
    return(NULL)
  }

  Npost <- 2000
  Nbin <- 2
  HaploidInd <- ifelse( haploid, 0, 1)

  if ( HaploidInd == 0 ) {
    model = "additive"
    Ndip = length(condensed.db$strains)
  }
  else {
    model = "full"
    ns = length(condensed.db$strains)
    Ndip = ns*(ns+1)/2
  }

  int = intersect( condensed.db$subjects, unique(phen.data$CC.Line))
  cat( length(int), "subjects analysed with phenotypes\n")

  phen.data = phen.data[phen.data$CC.Line %in% int,]
  db.use = match( phen.data$CC.Line, condensed.db$subjects)
  Nind = nrow(phen.data)
  
  result = condensed.db[[model]]$condensed.summary

  chr = condensed.db[[model]]$chr
  plps = mclapply(chr, scan.chromosome.hbrem, db.use, phen.data$x, HaploidInd, Ndip, Nind, Npost, Nbin, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores=mc.cores )
#  plps = lapply(chr, scan.chromosome.hbrem, db.use, phen.data$x, HaploidInd, Ndip, Nind, Npost, Nbin)

  res = t(do.call( "cbind", plps))

  mark.pars.df <- data.frame(res[,1:46])
  names(mark.pars.df) = c( "F.logPval", "Hbar", "sd.Ni", "BIC.qtl", "BIC.null", "BF", "logBF", "DIC.qtl", "DIC.null", "DIC.diff", "pd.qtl", "pd.null", "mode.kT", "ga", "gb", "mode.var", "med.kT", "med.mu", "med.var", "mean.kT", "mean.mu", "mean.var", "hpd.kT.50.lower", "hpd.kT.50.upper", "hpd.mu.50.lower", "hpd.mu.50.upper", "hpd.var.50.lower", "hpd.var.50.upper", "hpd.kT.75.lower", "hpd.kT.75.upper", "hpd.mu.75.lower", "hpd.mu.75.upper", "hpd.var.75.lower", "hpd.var.75.upper", "hpd.kT.95.lower", "hpd.kT.95.upper", "hpd.mu.95.lower", "hpd.mu.95.upper", "hpd.var.95.lower", "hpd.var.95.upper", "hpd.kT.99.lower", "hpd.kT.99.upper", "hpd.mu.99.lower", "hpd.mu.99.upper", "hpd.var.99.lower", "hpd.var.99.upper")

  offset = ncol(mark.pars.df)
  mark.pars.df$Name=condensed.db[[model]]$genome$marker 
  mark.pars.df$Chr = condensed.db[[model]]$genome$chromosome
  mark.pars.df$Bp = condensed.db[[model]]$genome$bp

  bp2 = mark.pars.df$Bp[2:length(mark.pars.df$Bp)]
  bidx = which(bp2 < mark.pars.df$Bp[1:(length(mark.pars.df$Bp)-1)])-1

  mark.pars.df$CumBp = rep(0,nrow(mark.pars.df))
  mark.pars.df$CumBp[bidx+2] = mark.pars.df$Bp[bidx+1]
  mark.pars.df$CumBp = cumsum(mark.pars.df$CumBp) + mark.pars.df$Bp
         
  strain.means.df <- data.frame( res[,offset:(offset+Ndip-1)])
  strain.sdevs.df <- data.frame( res[,(offset+Ndip):(offset+2*Ndip-1)])
  strain.avNis.df <- data.frame( res[,(offset+2*Ndip):(offset+3*Ndip-1)])

  if ( HaploidInd ==0 ) {
   names(strain.means.df) <- condensed.db$strains
   names(strain.sdevs.df) <- condensed.db$strains
   names(strain.avNis.df) <- condensed.db$strains
  }
  else {
   strain.names = condensed.db$strains
   num.strains = length(condensed.db$strains)
   diplotype.names <- matrix(kronecker(strain.names, strain.names, paste, sep="."), nrow=num.strains)
   names.full.symmetric <- c( diag(diplotype.names), diplotype.names[upper.tri(diplotype.names, diag=FALSE)])
   names(strain.means.df) <- names.full.symmetric
   names(strain.sdevs.df) <- names.full.symmetric
   names(strain.avNis.df) <- names.full.symmetric
  }


  hbrem.region.list <- list(Summary.Parameters=mark.pars.df, Group.Means=strain.means.df, Group.StDevs=strain.sdevs.df, Group.ExpCounts=strain.avNis.df)

  return(hbrem.region.list)

  }  	

scan.hbrem <- function( condensed.db, data.file, phen.names, haploid=FALSE, mc.cores=10, subject.name="CC.Line" ) {

  phen.data = read.delim(data.file)
  err=0
  for ( p in c( subject.name, "Sex", phen.names )) {
    if ( is.null( phen.data[[p]] ) ) {
      warning( "file ", data.file, " has no column ", p, "\n")
      err = err+1
    }
  }

  if ( err > 0)  return(NULL)

  scans = list()
  for ( phen.name in phen.names ) {
    if ( is.numeric(phen.data[[phen.name]]) ) {
     pdata.sex = data.frame( y=phen.data[[phen.name]], Sex=phen.data$Sex)
     cc = complete.cases(pdata.sex)
     pdata.sex = pdata.sex[cc,]
     agg.phen.data = aggregate( pdata.sex$y, by=list(CC.Line=phen.data[[subject.name]][cc], Sex=pdata.sex$Sex), FUN="mean")#FUN="median" )
     f = lm( agg.phen.data$x ~ agg.phen.data$Sex, data=agg.phen.data )
     agg.phen.data2 = aggregate( residuals(f), by=list(CC.Line=agg.phen.data$CC.Line), FUN="mean" )
     null.model = "x ~ 1"
     scan.lm = multicore.condensed.genome.scan.lm( condensed.db, agg.phen.data2, null.model, permute=1000, mc.cores=mc.cores)
     write.table(scan.lm$result, file=paste(phen.name, "lm.summary.txt", sep=""), quote=FALSE, row=FALSE, sep="\t")

     scan = multicore.condensed.genome.scan.hbrem( condensed.db, agg.phen.data2, haploid=haploid, mc.cores=mc.cores)
     write.table(scan$Summary.Parameters, file=paste(phen.name, ".hbrem.summary.txt", sep=""), quote=FALSE, row=FALSE, sep="\t")
     scans[[phen.name]] = scan
   }
  }
  return(scans)
}


plot.scans.hbrem <- function ( scans, pdffile=NULL ) {

  if ( ! is.null(pdffile)) pdf( pdffile ) 
  par(mfrow=c(4,1),cex=0.4)	

  for( p in names(scans)) {
    scan = scans[[p]]$Summary.Parameters
    chr = scan$Chr[2:nrow(scan)]
    boundary = c(1, which(chr != scan$Chr[1:nrow(scan)-1]))
    boundary2 = c(boundary[2:length(boundary)],nrow(scan))
    xc = 0.5*(boundary+boundary2)
    plot( scan$CumBp, scan$F.logPval, t="l", main=p, ylab="logP", xaxs="i", xaxt="n", xlab="position")
    abline( v=scan$CumBp[boundary], col="red")
    axis(1, at=scan$CumBp[xc], labels=chr[boundary])

    plot( scan$CumBp, scan$hpd.kT.50.upper, t="l", col="grey", main=p, ylab="mode.kT", xaxs="i", xaxt="n", xlab="position")
    lines( scan$CumBp, scan$hpd.kT.50.lower, t="l", col="grey")
    lines( scan$CumBp, scan$mode.kT, t="l", col="black")
    abline( v=scan$CumBp[boundary], col="red")
    axis(1, at=scan$CumBp[xc], labels=chr[boundary])
  }
  if ( ! is.null(pdffile )) dev.off()

}




# SURVIVAL ANALYSIS

scan.chromosome.survival <- function( mats, db.use, genetic.model, phen.data, fit.null, subjects ) {
  return( sapply( mats, test.locus.survival, db.use, genetic.model, phen.data, fit.null, subjects ))
}

test.locus.survival <- function( e, db.use, genetic.model, phen.data, fit.null, subjects ) {


  d = e$mat
  rownames(d) = subjects
  d = d[db.use,]
  phen.data$d = d
  fit.genetic = survreg( as.formula(genetic.model), data=phen.data)
  
  a = anova(fit.null, fit.genetic)
  
                                        #    return( c(-log10(a[2,6]), p))
  return( -log10(a[2,7]))
}

multicore.condensed.genome.scan.survival <- function( condensed.db, phen.data, null.model, permute=0, model="additive", mc.cores=8 ) {

  if ( class(condensed.db) != "condensed.happy" ) {
    warning( "condensed.db not of class condensed.happy\n")
    return(NULL)
  }


  int = intersect( condensed.db$subjects, unique(phen.data$CC.Line))
  cat( length(int), "lines with phenotypes and genotypes\n")
  cat( int, "\n")

  browser()
  phen.data = phen.data[phen.data$CC.Line %in% int,]
  line.model = as.formula( "Surv(phen.data$x,phen.data$death) ~ CC.Line + Sex")
  null.model = as.formula( "Surv(phen.data$x,phen.data$death) ~ Sex")
  base.model = as.formula( "Surv(phen.data$x,phen.data$death) ~ 1")
#  db.use = phen.data$CC.Line
  db.use = match( phen.data$CC.Line, condensed.db$subjects)


  fit.line = survreg( line.model, data=phen.data );
  fit.null = survreg( null.model, data=phen.data );
  fit.base = survreg( base.model, data=phen.data );
  print(anova(fit.base,fit.null,fit.line,test="Chisq"))

  genetic.model = as.formula( "Surv(phen.data$x,phen.data$death) ~ Sex + d")
  result = condensed.db[[model]]$condensed.summary

  chr = condensed.db[[model]]$chr
#  plps = lapply(chr, scan.chromosome.survival, db.use, genetic.model, phen.data, fit.null, condensed.db$subjects)
  plps = mclapply(chr, scan.chromosome.survival, db.use, genetic.model, phen.data, fit.null, condensed.db$subjects, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores=mc.cores )

  logP = NULL
  for ( i in 1:length(plps))
      logP = c( logP, plps[[i]])
  
  result$logP = logP
  max.logP = max(logP)
  cat("genome scan max logP", max.logP, "\n")
  perm.logP.sorted = NULL
  cat("permute ", permute, "\n")
  if ( permute > 0 ) {
	perm.logP = sapply( 1:permute, function( p, condensed.db, db.use, genetic.model, phen.data, fit.null, subjects, max.logP, mc.cores=mc.cores  ) {

#	  db.use.perm = sample(db.use,replace=FALSE)
  	  db.use.perm = match( phen.data$CC.Line, sample(condensed.db$subjects,replace=FALSE))
	  chr = condensed.db[[model]]$chr
	  plps.perm = mclapply(chr, scan.chromosome.survival, db.use.perm, genetic.model, phen.data, fit.null, subjects, mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores=mc.cores )

	  perm.logP = NULL
	  for ( i in 1:length(plps))
      		perm.logP = c( perm.logP, plps.perm[[i]])
	  perm.logP.max = max(perm.logP)
	  cat(p, "genome scan permuted max logP", perm.logP.max, "\n")
          return(perm.logP.max)
	}, condensed.db, db.use, genetic.model, phen.data, fit.null, condensed.db$subjects, max.logP, mc.cores=mc.cores  )
	perm.logP.sorted = sort(perm.logP)
        
	}

	

  return( list(result=result, permuted.logP=perm.logP.sorted ) )
}  	

survival.median <- function(y) { 
  quantile( y, probs=0.5, type=1, na.rm=TRUE); 
}


scan.survival <- function( condensed.db, data.file, surv.time, surv.censored=28, permute=0, histogram.pdf=NULL, mc.cores=10, offset=0 ) {

  cat("offset=",offset,"surv.censored=",surv.censored,"\n")
  phen.data = read.delim(data.file)
  null.model = "~ Sex"
  err=0
  for ( p in c( "CC.Line", "Sex", surv.time )) {
    if ( is.null( phen.data[[p]] ) ) {
      warning( "file ", data.file, " has no column ", p, "\n")
      err = err+1
    }
  }

  if ( err > 0)  return(NULL)

  scans = list()
  if ( is.numeric(phen.data[[surv.time]]) ) {

   su = Surv(phen.data[[surv.time]]+offset,phen.data[[surv.time]]<surv.censored+offset)
   complete.null = survreg( su ~ 1 )
   complete.sex = survreg( su ~ Sex, data=phen.data)
   complete.sex.line = survreg( su ~ Sex + CC.Line, data=phen.data )
   print(anova(complete.null,complete.sex,complete.sex.line))

   agg.phen.data = aggregate( phen.data[[surv.time]], by=list(CC.Line=phen.data$CC.Line, Sex=phen.data$Sex), FUN="mean")#FUN="median" )
   agg.phen.data$death = agg.phen.data$x < surv.censored
   s = Surv( agg.phen.data$x, agg.phen.data$death )
#  scan = condensed.genome.scan.lm( condensed.db, agg.phen.data, paste( , null.model), permute=permute, phen.name=phen.name)
   scan = multicore.condensed.genome.scan.survival( condensed.db, agg.phen.data,  null.model, permute=permute, mc.cores=mc.cores)
   write.table(scan$result, file=paste(p, ".scan.txt", sep=""), quote=FALSE, row=FALSE, sep="\t")
     scans[[p]] = scan
  }
  phen.data.plot = data.frame(x = agg.phen.data$x)
  names(phen.data.plot) = c( surv.time)
  plot.histograms( histogram.pdf, phen.data.plot, c(surv.time) ) 

  return(scans)
}

cousin.lines <- function( ids, rm=TRUE ) {
	splat = strsplit( as.character(ids), "[-.]", extended=TRUE )
	d = data.frame(do.call( "rbind", splat ) )
	names(d) = c("prefix", "N")
	d$id = ids;
	d$N = as.numeric(as.character(d$N))
	d$line = ifelse( d$N> 500, d$N-500 , d$N)
	
	r = order(d$line)
	df = d[r,]
	if ( rm==TRUE ) {
		u = unique(df$line)
		df = df[match(u, df$line ),]
	}
	return(df)
}

plot.haplotypes.in.region <- function( g, subjects, chr, from, to, file="haplotypes.pdf" ) {

	pdf(file)
	s = match( subjects, g$subjects, nomatch=0)
	ns =length(subjects)
	n = max(4,ns) 
	par(mfrow=c(4,1))
	w = g$additive$genome$chromosome==chr & g$additive$genome$bp >= from & g$additive$genome$bp <= to
	markers = g$additive$genome$marker[w]
	nm = length(markers)
	
	mat = list()
	for ( sub in subjects ) 
		mat[[sub]] = matrix(0,nrow=8,ncol=nm)
	i = 1
	for( m in markers ) {
		d=condensed.hdesign(g, m)
                ds = apply(d,1,sum)
                d = d/ds
		for( j in 1:ns )
			mat[[subjects[j]]][,i] = d[s[j],]
		i = i+1
	}
	
	labels = g$strains
	for ( sub in subjects ) {
		image( t(mat[[sub]]), axes=FALSE, col=grey(1-(1:10/10)), main=sub)
		axis(2,at=seq(0,length(labels)-1)/(length(labels)-1),labels=labels,las=1,cex.axis=0.6)
		axis(1,at=seq(1,nm)/nm,labels=sprintf("%.2f", g$additive$genome$bp[w]/1.0e6))
	}

	dev.off()
}

plot.haplotypes <- function( g, subjects=g$subjects, chrs=paste( "chr", c(1:19, "X"), sep=""), file="haplotypes.pdf", thin=5, subject.file=NULL ) {

    texts = NULL
    subjects.txt = subjects
    if ( !is.null(subject.file)) {
        texts = read.table(subject.file, header=TRUE)
        subjects = texts$CC.Line
        subjects.txt = paste( texts$type, texts$CC.Line )
    }
   
    pdf(file)
    subjects = intersect(subjects,g$subjects)
    s = match( subjects, g$subjects, nomatch=0)
    ns =length(subjects)
    n = max(4,ns) 
    par(mfrow=c(4,1))
    w = g$additive$genome$chromosome %in% chrs
    cat( "nm: ", length(w), "\n")
    w = w[seq(1,length(w),thin)]
    cat( "nm.thin: ", length(w), "\n")
    markers = g$additive$genome$marker[w]
    nm = length(markers)
    bp = g$additive$genome$bp[w]
    bound = which( diff(bp) < 0 )
    print(bound)
    mat = list()
    for ( sub in subjects ) 
        mat[[sub]] = matrix(0,nrow=8,ncol=nm)
    i = 1
    for( m in markers ) {
        d=condensed.hdesign(g, m)
        for( j in 1:ns )
            mat[[subjects[j]]][,i] = d[s[j],]
        i = i+1
    }
    
    labels = g$strains
    k = 0
    for ( sub in subjects ) {
        k = k+1

        image( t(mat[[sub]]), axes=FALSE, col=grey(1-(1:10/10)), main=subjects.txt[k])
        
        axis(2,at=seq(0,length(labels)-1)/(length(labels)-1),labels=labels,las=1,cex.axis=0.6)
                                        #		axis(1,at=seq(1,nm)/nm,labels=sprintf("%.2f", g$additive$genome$bp[w]/1.0e6))
        abline( v=bound/nm, col="red")
    }
    
	dev.off()
}	


		
	
