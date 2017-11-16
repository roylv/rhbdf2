library(plotrix)

bone.h2 <- function( file="bone.v1.csv", cols=c("CORT.Dia.Dia", "CORT.Med.Dia", "CORT.Cort.Th", "CORT.Cort.Por", "CORT.BS.BV", "CORT.Mean1", "MOI.Polar", "MOI.Areal", "FULL.BV.TV", "FULL.Conn.D", "FULL.SMI", "FULL.Tb.N", "FULL.Tb.Th", "FULL.Sp", "PROX.BV.TV", "PROX.Conn.D", "PROX.SMI", "PROX.Tb.N", "PROX.Tb.Th", "PROX.Sp", "DIST.BV.TV", "DIST.Conn.D", "DIST.SMI", "DIST.Tb.N", "DIST.Tb.Th", "DIST.Sp", "Imax.Imin", "Length"),pdf.file="bone.hist.pdf")  {

  b = read.csv(file)
  pdf(pdf.file)
  par(mfrow=c(3,3))
  for( x in cols ) {
    cc = complete.cases( b$SUBJECT.NAME, b$Sex, b[[x]])

    f00 = lm ( b[[x]] ~ 1, data=b, subset=cc)
    f0 = lm ( b[[x]] ~ b$Sex, data=b, subset=cc)
    f1 = lm ( b[[x]] ~ b$Sex + b$SUBJECT.NAME, data=b, subset=cc)
    a = anova( f00, f0, f1)
    h = (a[2,2]-a[3,2])/(a[2,2])
    h.sex = (a[2,2]-a[3,2])/(a[1,2])
    sex = (a[1,2]-a[3,2])/(a[1,2])
    cat( x, h, a[3,6], h.sex,
        sex, a[2,6],  "\n")
    print(a)
    male = b[[x]][cc & b$Sex =="M"]
    female =b[[x]][cc & b$Sex =="F"]
    multhist( list( male=male, female=female), main=x)
    cat("\n\n")
  }
  dev.off()
}

      
    
  
