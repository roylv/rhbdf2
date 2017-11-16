nomda=read.delim("D:/analysis70/mf/7jan15/trab/BVTV.scanmarch15.txt")
mda=read.delim("D:/analysis70/mf/7jan15/trab/BVTV_MDA.scan.txt")
mda2=read.delim("D:/analysis70/mf/7jan15/trab/BVTV_MDA2.scan.txt")
nomda4=read.delim("D:/analysis70/mf/7jan15/trab/BVTV_4.scan.txt")
nomda5=read.delim("D:/analysis70/mf/7jan15/trab/BVTV_5.scan.txt")

par(mfrow=c(1,1))
plot(mda$logP,type="l",lwd=0.5)
lines(nomda$logP,type="l",lwd=0.5,col="red")
chr2mda=mda[mda$chromosome=="chr2" & mda$from.bp>118576459 & mda$from.bp<132587523,]
chr2mda2=mda2[mda2$chromosome=="chr2" & mda2$from.bp>118576459 & mda2$from.bp<132587523,]
chr2nomda4=nomda4[nomda4$chromosome=="chr2" & nomda4$from.bp>118576459 & nomda4$from.bp<132587523,]
chr2nomda5=nomda5[nomda5$chromosome=="chr2" & nomda5$from.bp>118576459 & nomda5$from.bp<132587523,]

str(chr2mda)
plot(chr2nomda$from.bp,chr2nomda$logP,type="l",ylim=c(0,8),xlim=(c(1.18e+08,1.32e+08)))
chr2nomda=nomda[nomda$chromosome=="chr2" & nomda$from.bp>118576459 & nomda$from.bp<132587523,]
lines(chr2mda$from.bp, chr2mda$logP,type="l",col="blue",lwd=2.5)
lines(chr2mda$logP,type="l",col="red",lwd=2.5)
lines(chr2mda2$from.bp, chr2mda2$logP,type="b",col="purple",lwd=2)
quantiles=c( 3.043246 ,4.028873, 4.273416)


lines(chr2nomda$from.bp,chr2nomda$logP,type="p",ylim=c(0,8),xlim=(c(1.18e+08,1.32e+08)))

pdf("noMDAvsMDAchr2_3.pdf", width=9.19, height=6.83)

plot(chr2nomda$from.bp,chr2nomda$logP,type="b",ylim=c(0,8),xlim=(c(1.18e+08,1.32e+08)),lwd=1.5)
lines(chr2nomda$from.bp,chr2nomda$logP,type="l",ylim=c(0,8),xlim=(c(1.18e+08,1.32e+08)),lwd=1.5)
lines(chr2mda$from.bp, chr2mda$logP,type="b", col="blue",lwd=1.5)
lines(chr2mda$from.bp, chr2mda$logP,type="l", col="blue",lwd=1.5)
lines(chr2mda2$from.bp, chr2mda2$logP,type="b",col="purple",lwd=1.5)
lines(chr2mda2$from.bp, chr2mda2$logP,type="l",col="purple",lwd=1.5)

lines(chr2nomda4$from.bp, chr2nomda4$logP,type="b",col="coral3",lwd=1.5)
lines(chr2nomda4$from.bp, chr2nomda4$logP,type="l",col="coral3",lwd=1.5)
lines(chr2nomda5$from.bp, chr2nomda5$logP,type="b",col="green",lwd=1.5)
lines(chr2nomda5$from.bp, chr2nomda5$logP,type="l",col="green",lwd=1.5)

legend("topright",c("No MDA","MDA","MDA no 2438","MDA only 2438","No MDA No 2438"),lty=c(1,1,1,1,1),lwd=c(1,2.5,2.5,2.5,2.5),col=c("black","blue","purple","coral3","green"))
abline(h=quantiles,col="red")
dev.off()
getwd()



#merge
mrg=read.delim("D:/analysis70/mf/7jan15/trab/BVTV.chr1_26700000-34200000.merge.txt",sep="")
names(mrg)
head(mrg$entropy)
plot(mrg$entropy,mrg$logP.merge)

a=mean(mrg$entropy)-sd(mrg$entropy)
b=mean(mrg$entropy)+sd(mrg$entropy)

plot(mrg$entropy,exp(-mrg$logP.merge),col=ifelse(mrg$entropy<a, "red", ifelse(mrg$entropy>b, "blue","black")),lwd=2,pch=ifelse(mrg$entropy<a, 1, ifelse(mrg$entropy>b, 2,3)),xlim=c(-1.5,-0.5))
abline(v=mean(mrg$entropy),col="coral",lwd=2)
abline(v=median(mrg$entropy),col="coral3",lwd=2)
abline(h=mean(mrg$logP.merge),col="coral3",lwd=2)

mean(mrg$entropy)


trab=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/PrelimTrab.txt")
names(trab)
nonvarst=c(1,2,4,5,6,7,8)
trab2=trab[,-nonvarst]
names(trab2)
trabLog=glm(Sex~.,data=trab2,family=binomial)

trabLog2=step(trabLog)
summary(trabLog2)  
Sex2=sapply(trab2$Sex,function(x) if(x=="M") x="1" else x="0")
trab$Sex2=Sex2
trab2=trab2[,-1]
str(trab2)
names(trab)
trab2$Sex=Sex2

trab2$Sex=as.numeric(trab2$Sex)
levels(trab2$Sex)=levels(factor(trab2$Sex))
trabFit=lm(Length~.,data=na.omit(trab2))
summary(trabFit)
plot(trabFit)
cor(trab2$Length,trab2$SMID)
trabFit2=step(trabFit)

summary(trabFit2)

library(rafalib)
mypar(1,1)
plot(trab2$Length,trab2$Med.Dia,type="p")
names(coef(trabFit2))
fitT=lm(Med.Dia~Length ,data=trab2)

abline(coef(fitT)[1],coef(fitT)[2],col="black")
coef(fitT)
predict(fitT, new)

coef(fitT)["Med.Dia"]

fitT2=lm(Mean2~Length ,data=trab2)
plot(trab2$Length,trab2$Mean2,type="p")
abline(coef(fitT2)[1],coef(fitT2)[2],col="black")

plot(trabLog2)

summary(trab2$Med.Dia)

full=read.csv("D:/analysis70/mf/7jan15/trab/fullmarch15.csv")
full2=read.csv("D:/analysis70/mf/7jan15/trab/allmarch152.txt",header=T,sep="\t")
summary(full2)
full3=full2
names(full2)
#fill-in NAs
library(mice)
set.seed(144)

unique(rownames(full2))
summary(loan)
vars.for.imputation = setdiff(names(full2), names(full2)[c(-33,-18)])
unique(vars.for.imputation)
imputed = complete(mice(full2[vars.for.imputation]))
full2[vars.for.imputation] = imputed

full3[,"weight"]
imputed = complete(mice(full2[vars.for.imputation]))
full2[vars.for.imputation] = imputed
ncol(full)
for(i in 1:nrow(full2)) {for(j in 1:nrow(full2)) if (full2[j,i]=="-") full2[j,i]=NA}

for(i in 1:nrow(full2)) {if(full2[i,"weight"]==0) full2[i,"weight"]=NA}
nrow(full2)

names(full2)

mean(!is.na(full2[,"weight"]))

colnames(full2)[15]

fit=lm(full3$Length~full3$BV.TV.met)

coefs=coef(fit)
coef(fit)
coefs
plot(full3$BV.TV.met,full3$Length)

abline(coefs[1],coefs[2])

names(full2[13:40])

rnorm(1)


getwd()
write.csv(full2,file="D:/analysis70/mf/7jan15/trab/march15All.csv")
?write.csv

(full2$BV.TV.met)
t.test( c(0.14,0.23),full2$BV.TV.met[1:nrow(full2)])

mean(c(0.13,0.2))



tb=read.csv("D:/analysis70/mf/7jan15/trab/march15All.csv")
summary(tb)
any(is.na(tb))

na=sapply(1:ncol(tb), function(x) if (any(is.na(tb[,x]))) return(colnames(tb)[x]))

na=sapply(1:ncol(tb), function(x) any(is.na(tb[,x])))

names(tb)[-which(na)]

colnames(tb[na])
which(na)
length(na)

library(mice)

set.seed(144)

unique(rownames(full2))
summary(loan)
vars.for.imputation = setdiff(names(tb), names(tb)[-which(na)])
unique(vars.for.imputation)
imputed = complete(mice(tb[vars.for.imputation]))
tb[vars.for.imputation] = imputed

write.csv(tb,file="D:/analysis70/mf/7jan15/trab/imputed15.csv")
names(tb)

library(rafalib)

names(tb)

setwd("C:/Users/user/Google Drive/YG/CC Proj/redHawk")

pdf("qqplots_march15.pdf",width=8.27, height=11.69)
mypar2(6,6)
for(i in c(16:ncol(tb), 12))  {qqnorm(as.numeric(tb[,i]),main=names(tb)[i]); qqline(as.numeric(tb[,i]))}

dev.off()

sapply(1:ncol(tb), function(x) .factor(tb[,x]))

tb$BV.TV.met=as.factor(tb$BV.TV.met)

qqnorm(as.numeric(tb$BV.TV.met))

qqnorm(tb$Areal.cort)
49-16

pdf("boxplots_march15.pdf",width=8.27, height=11.69)
mypar2(6,6)
for(i in c(16:ncol(tb), 12))  {boxplot(as.numeric(tb[,i]),main=names(tb)[i])}

dev.off()


spec=tb$SPECIMEN
spec
spec=sub("(.*?)_FEMUR.*", "\\1", spec)
spec=sub("(.*?)_femur.*", "\\1", spec)
spec=gsub("CC#","", spec)
spec=gsub("cc#","", spec)
spec=gsub("IL-","", spec)

spec=tb$CC.Line

spec=as.numeric(spec)

spec
spec=paste("IL-",spec,sep="")

tb$CC.Line=spec
tb$CC.Line

write.csv(tb,file="D:/analysis70/mf/7jan15/trab/imputed15.csv")
names(tb)
tb$line.Label
tb2=tb
spec2=spec

tb2$CC.Line=spec2

for(i in 1:nrow(tb2)) if(tb2[i,"line.Label"]=="B") {tb2[i,"CC.Line"]=500+tb2[i,"CC.Line"]}

tb2$CC.Line=paste("IL-",spec,sep="")

tb2$CC.Line=as.factor(tb2$CC.Line)
tb$CC.Line=tb2$CC.Line

tb$CC.Line
a=84*2+98*2+91*4+83*9+78+98+95*2+93
a/22
names(tb[c(16:49,12)])

write.table(tb,file="imputed15.txt",sep="\t",quote=FALSE,row=FALSE)
getwd()


fc <- file("http://mus.well.ox.ac.uk/preCC/MEGA_MUGA/Mar2015.MEGA+MDA+MUGA/chr1.MEGA+MDA+MUGA.alleles")
mylist <- strsplit(readLines(fc), "\t")
close(fc)
data.frame(mylist)


#alleles
getwd()

write.table(rd3,file="c:/temp/allele4.txt",sep="\t",quote=FALSE,row=FALSE)
rdl=read.delim("c:/temp/allele4.txt",header=T)
rdl
class(rd3)
rdl=read.delim(file="clipboard",header=T)
rd3=na.omit(rdl)
rd3=rd3[,-1]

1194%%4
for(i in 1:nrow(rd2)) if(as.numeric(rownames(rd2[i,]))%%2==1) rd2[i,]=rd2[i,]*-1 
rd2
seq=(1:nrow(rdl))
rdl=rdl[-which(seq%%4==0),]
sm2=apply(rd2,2,sum)
sm2-sm2["V9"]
barplot(sm2-sm2["V9"])/sum(sm2)





#HEATMAP EXPERIMENT
library(matrixStats)
dat=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/imputed/data28.txt")
info=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/imputed/info28.txt")
group=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/imputed/group28.txt")
dat=as.matrix(dat)
group=factor((as.list(group))$x)

names(info)
top50=order(-rowMads(dat))[1:25]

top50=order(rowttests(dat),decreasing=T)[1:25]
top50=order(apply((t(dat)),2,sd),decreasing=T)[1:35]


nm=c('bv1f','bv2f','bv3f','tb1f','tb2f', 'tb3f','conn1f','conn2f','conn3f'
     ,'bv1m','bv2m','bv3m','tb1m','tb2m', 'tb3m','conn1m','conn2m','conn3m')

phenGroup=c(1,1,1,1,1,1,3,3,3,3,3,3,1,1,1,3,3,3)
phenGroup=factor(phenGroup)
colnames(dat)

library(rafalib)
?heatmap.2
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
cols=palette(brewer.pal(11,"RdBu"))[phenGroup]
heatmap.2(dat[top50,],
          #labCol=sampleInfo$date,
          trace="none"
          ,ColSideColors=cols,
          col=hmcol,
          labRow=info$chromosome[top50],
          scale="row",
          key=F)
 head(sort(dat[,18],decreasing=T))
top50





##boxplots
setwd("YG/CC Proj/analyses")
dat=read.csv("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/march15/march15ALL.csv")
names(dat)
dat2=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/march15/imputed15.txt")
dat$CC.Line=dat2$CC.Line

sx=split(dat, dat$Sex)
bvm=split(sx$M$BV.TV.met, sx$M$CC.Line)
bvf=split(sx$F$BV.TV.met, sx$F$CC.Line)
mypar(2,1)
?boxplot
pdf("boxplotsBVTVALL.pdf")
mypar(3,1)
boxplot(bvm, col="blue", range=0, cex.axis=0.5,las=2)
boxplot(bvf, col="pink", range=0, cex.axis=0.5,las=2)
dev.off()

install.packages("outliers")
library(outliers)
outF=outlier(sx$F$BV.TV.met,logical=T)

sx$F=sx$F[-which(outF),]

remove_outliers <- function(x, na.rm = TRUE) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
sx$M$BV.TV.met=remove_outliers(sx$M$BV.TV.met)

qnorm(sx$M$BV.TV.met)
qnorm

datSort=dat[order(dat$BV.TV.met),]

plot(datSort$BV.TV.met)


boxplot(datSort$BV.TV.met~datSort$CC.Line, range=0.97, las=2, cex.axis=.5, col="darkgray")
boxplot(datSort$BV.TV.dist~datSort$CC.Line, range=0.97, las=2, cex.axis=.5, col="darkgray")
mypar(3,1)

(bvm, col="blue", cex.axis=0.5,las=2)


df <- data.frame(f1=factor(rbinom(100, 1, 0.45), label=c("m","w")), 
                 f2=factor(rbinom(100, 1, 0.45), label=c("young","old")),
                 boxthis=rnorm(100))

df$f1f2 <- interaction(df$f1, df$f2)
df$f2
dat$CC.LineMF=interaction(dat$Sex, dat$CC.Line)

gp=ggplot(aes(y = BV.TV.met, x = CC.LineMF, fill=Sex), data = dat) + geom_boxplot() + scale_fill_manual(values=c("blue","pink"))
gp+ theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.background = element_rect(fill = "white",colour = NA), # or theme_blank()
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "white",colour = NA)
  )


data(diamonds)
names(diamonds)

df$boxthis
mypar(1,1)
library("vioplot")
boxplot(BV.TV.met~Sex*CC.Line, data=dat, notch=F, col=(c("pink", "blue")), las=2, range=0.5, cex.axis=.5)

barplot(sx$F$BV.TV.met, col="lightblue")
lines(density(sx$F$BV.TV.met))


ft=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/MCRemoved/Trab/MCSim.txt")
dim(ft)
ftbv=ft[which(ft[,11]=="BV.TV.met"),]
ftth=ft[which(ft[,11]=="Tb.Th.met"),]
fttb=ft[which(ft[,11]=="Tb.N.met"),]
dim(ftbv)
plot(fttb[fttb[,1]>0,1],type="l", axes=F, ylab="", xlab="",pch=19); axis(side=4)
lines(fttb[,4], col="red")
par(new=T)   #tell R to overwrite the first plot
plot(fttb[,6], col="blue", type="p", pch=19, ylab="", xlab="")
abline(v=91)

fttb[which.max(fttb[,1]),1]

levels(factor(dat$D.B))
hist((dat$BV.TV.met))
lines

sum(dat$Generat>20)
dte=read.table(file="clipboard")
db=(as.list(dat$D.B2))
dat$D.B2=as.Date(dat$D.B2)

getwd()
write.table(dat,file="march15ALLdates.txt",sep="\t",quote=FALSE,row=FALSE)
dat=read.delim("YG/CC Proj/redHawk/march15ALLdates.txt")
class(dat$D.B2)
dat2$DB=dat$D.B2
dat2=dat
names(dat)
yr=format(as.Date(dat$D.B2), "%y")
mnt=format(as.Date(dat$D.B2), "%m")
hist(as.numeric(mnt))
levels(factor(yr))
which(mnt=="06" | mnt=="07")
dat2=dat2[-which(mnt=="06" | mnt=="07"),]

length(which(table(dat2$CC.Line)>0))

names(dat2)

table(mnt)
dat2$DB=as.Date(dat2$DB)
mnt2=format(dat2$DB, "%m")
table(factor(mnt2))
qqnorm(dat2$Mean1.met)

summary(lm(dat$Cort.Th~mnt, data=dat))
factor(dat$Week.old)

mnt


names(dat)
month=format(as.Date(dat$D.B2), "%m")
monthyear=format(as.Date(dat$D.B2), "%m%y")
year=format(as.Date(dat$D.B2), "%y")

month

names(dat)

dat$Tb.N.prox=as.numeric(dat$Tb.N.prox)
boxplot(split(dat$Tb.N.met, day),col="darkgreen",las=2)

mont=format(as.Date(dat$D.B2), "%m")
day=format(as.Date(dat$D.B2), "%d")
  
idxyr=which(yr=="13")
datyr=dat[idxyr,]
mntyr=month[idxyr]

sp=split(dat$Tb.N.met, daymnt)

plot(sapply(sp, function(x) mean(x)), pch=19)

stripchart(split(dat$Tb.N.met, day),add=TRUE,vertical=TRUE,pch=19,cex=.5,col=1)
}
X

summary(lm(dat$BV.TV.met~X[,2]))

mean(table(month))

X=model.matrix(~mnt)

mnt=X[,c(2:3)]%*%matrix(c(1,1),2,1)
mnt

boxplot(split(dat$Tb.N.met, month),col="darkgreen", las=2, range=0.5)
mnsplt=split(dat$Tb.N.met, mnt)
a=mean(mnsplt[[1]])
b=mean(na.omit(mnsplt[[2]]))
barplot(a,b)


sapply(mnsplt, function(x) mean())
barplot()
library(genefilter)
idx1=dat$BV.TV.met[mnt==1]
idx2=dat$BV.TV.met[mnt==0]

t.test(dat$BV.TV.met[idx1],dat$BV.TV.met[idx2])

library(genefilter)
rowttests(t(matrix(dat$Tb.Th.met[-idx])),factor(mnt[-idx]))
month

idx=which(is.na(dat$Tb.Th.met))

dat$mnt=mnt
fit=lm(dat$BV.TV.met~dat$Sex+mnt)

table(dat$CC.Line, mnt)


idx2=which(dat$CC.Line %in% c("IL-2750","IL-3912","IL-4052"))

dat$Sex
plot(resid(fit))

summary(fit)
plot(dat$BV.TV.met, dat$weight)

cor(dat$BV.TV.met,mnt)

factor(month)
boneCART=rpart(BV.TV.met~, data=dat, method="class")

prp(boneCART)

randomForest(dat$BV.TV.met~mnt+dat$Sex+dat$Week.old)

setwd("YG/CC Proj/redHawk")
dat2=read.delim("imputed15.txt")
dat$CC.Line=dat2$CC.Line
fit0=lm(dat$BV.TV.met~1)
fit1=lm(dat$BV.TV.met~mnt)
fit2=lm(dat$Tb.N.met~mnt)
fit3=lm(dat$BV.TV.met~yearmnt)

summary(lm(dat$BV.TV.met~yearmnt))
a=anova(fit0,fit1)


mean(dat$BV.TV.met[which(mnt==0)])

(408-362)/408

h = (a[2,2]-a[3,2])/(a[2,2])

h.sex = (a[2,2]-a[3,2])/(a[1,2])

(0.8-0.7)/0.8

(0.76-0.33)/0.86
(dat$BV.TV.met[dat$CC.Line=="IL-557"])
mean(dat$Generat)

mnt

dat14=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/BV.TV.met.scan.txt")
dat146=dat14[dat14$chromosome=="chr6",]

dat146$from.bp[which.max(dat146$logP)]

sim=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/MCSim.txt")
simsub=sim[sim$V1>=6 & sim$V8>=39,]
simsub


names(table(dat$CC.Line))

plot(t(cor(as.numeric(as.factor(dat$weight)),dat[,c(20:31,33:43,45:53)])))

plot(dat$Length, dat$weight)

#pheno103
dat103=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/pheno103date.txt")
str(dat103)
nrow(dat103)
sum(complete.cases(dat103))
anova(lm(dat103$Tb.N3D~dat103$Sex))
names(dat103)
dat103$BMC
conn103=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/Conn.D.scan.txt")
connd14=conn103[(conn103$chromosome=="chr14"),]

connd14$to.bp[which.max(connd14$logP)]

#mcsim
mcsim=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/MCSim.txt")
which.max(mcsim$V1)
mcsim[3337,]
ln3337=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/linesRemoved3337.txt")
ln3337


#cort
ct=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/Cort.Th.scan.txt")
ct11=ct[which(ct$chromosome=="chr11"),]

ct11$from.bp[which.max(ct11$logP)]

pm=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/dateStrat/Cort.Th.scan.txt")
pm11=pm[which(pm$chromosome=="chr11"),]

pm11$to.bp[which.max(pm11$logP)]


dat=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/Tb.NPl.scan.txt")
dat1=dat[which(dat$chromosome=="chr1"),]

dat$to.bp[which.max(dat$logP)]


###chip
f2=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/ENCFF000KSW.bed", header=F)
f2i=import(f2, format="bedGraph")
class(f2)

which(f2$id=="OXT")
f2[21003,]
which.min(f2$score)
head(f2)

colnames(f2)=c('chr','start','end','id','score','strand')
f2b= with(f2, GRanges(chr, IRanges(start, end), strand, score, id=id))
summary(f2b)
class(f2b)
plotRanges(ranges(f2b))



pheno="Tb.NPl"
ylb="rand"
####barplot
library(plotrix)  
library(ggplot2)
descBars = function(df, ylb, fill="dark blue", pheno, desc=T){                                                                 

  df[[pheno]]=as.numeric(df[[pheno]])
  spdf2=split(df, df$CC.Line)
  se=sapply(spdf2, function(x) std.error(x[[pheno]]))
  mean(df[[pheno]])
  df[[pheno]]=as.numeric(df[[pheno]])
  mbvtvm=sapply(spdf2, function(x) mean(x[[pheno]]))
  
  
  fit=lm(df[[pheno]]~CC.Line, data=df)
  #nm=names(which(coef(fit)<0.001))
  #nm=gsub("CC.Line", "",nm)
  
  smfit=summary(fit)
  pv=smfit[4][[1]][,4]
  pv=pv[-1]
  nm1=gsub("CC.Line", "",names(which(pv<=0.05 & pv>0.01)))
  nm2=gsub("CC.Line", "",names(which(pv<=0.01)))
    
  #mbvtvm[nm1]
  #m=rep(0, length(nm))
  lb.df1=data.frame(l=nm1, m=mbvtvm[nm1]+1.1*se[nm1])
  lb.df2=data.frame(l=nm2, m=mbvtvm[nm2]+1.1*se[nm2])
  lb.df1.1 = data.frame(l=names(which(is.na(se[nm1]))), m=1.01*mbvtvm[names(which(is.na(se[nm1])))])
  lb.df2.1 = data.frame(l=names(which(is.na(se[nm2]))), m=1.01*mbvtvm[names(which(is.na(se[nm2])))])
  
  a=list()
  a$se=se
  a$se=a$se[order(names(a$se))]
  a$m=mbvtvm
  a$m=a$m[order(names(a$m))]
  a$l=sort(names(table(df$CC.Line)))

  a = as.data.frame(a)
  if(desc) {a = transform(a, l=reorder(l, -m))}
  p=ggplot(a, aes(x=l, y=m))+
    geom_bar(stat="identity", fill=fill) + 
    ylab(ylb) +
    theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,hjust=1))
  
  # Define the top and bottom of the errorbars
  limits <- aes(ymax = m + se, ymin=m)
  #p + geom_bar(fill=fill, position=position_dodge(), stat="identity") + geom_errorbar(limits)#, position=dodge, width=0.25)
  if(sum(is.na(lb.df1))==nrow(lb.df1) & sum(is.na(lb.df2))==nrow(lb.df2)){
    
    p+p+geom_errorbar(limits) + 
      geom_text(data = lb.df1.1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df2.1, label = "**", color="red", size=7)
  }
  else if(sum(is.na(lb.df1))==nrow(lb.df1) & sum(is.na(lb.df2))<nrow(lb.df2)){
    p+geom_errorbar(limits) + 
      geom_text(data = lb.df2, label = "**", color="red", size=7)+
      geom_text(data = lb.df1.1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df2.1, label = "**", color="red", size=7)
    
  }
  else if(sum(is.na(lb.df1))<nrow(lb.df1) & sum(is.na(lb.df2))==nrow(lb.df2)){
    p+geom_errorbar(limits) + 
      geom_text(data = lb.df1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df1.1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df2.1, label = "**", color="red", size=7)
    
  }
  else{
  
    p+geom_errorbar(limits) + 
      geom_text(data = lb.df1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df2, label = "**", color="red", size=7)+
      geom_text(data = lb.df1.1, label = "*", color="dark red", size=7)+
      geom_text(data = lb.df2.1, label = "**", color="red", size=7)
  }
  
}
library(ggplot2)
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/prelim103.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/FemOctNoOut.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/biom/BiomClean.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/Cal/CalComb2.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/Vert/Vert100.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/maxScan/Dec/BVTV3/BVTV258/BVTV258.txt")

df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/maxScan/Dec/BVTV3/BVTV258/BVTV258.txt")
df1=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/FemOctNoOut.txt")
setwd("c:/users/Roy/google drive/YG/CC Proj/rhbdf/")
load("matA_BV.TV.RData")
load("matR_BV.TV.RData")
ccla=mA[36675,][[1]]
cclr=mR[36675,][[1]]
df=rbind(df,df1[which(df1$CC.Line %in% ccla),])
df=read.delim("BVTV36675.txt")

dfm=df[which(df$Sex=="M"),]
dff=df[which(df$Sex=="F"),]

length(unique(df$CC.Line))

dfn=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/maxScan/Dec/BVTV3/BVTV258/BVTV436_2U.txt")



dim(df)[[1]]+29

length(which(df$CC.Line=="IL-557"))

df=df[-which(df$CC.Line %in% cclr),]

write.table(df,file="BVTV36675.txt",sep="\t",quote=FALSE,row=FALSE)


ccla %in% cclr

names(df)

phen=c(colnames(df)[c(7, 8, 10, 11)])
phen=c("vBMD","Porosity","BMC","Dia.Dia","Med.Dia","Cort.Th","Areal","BV.TV","Conn.D","SMI","Tb.N3D","Tb.NPl","Tb.Th", "Sp", "Mean1")
phen=c("BV.TV","Tb.NPl","Tb.Th","Conn.D","SMI","Sp", "BV/TV", "Tb.N", "Tb.Th", "Conn.D", "SMI", "Sp")
phen=c("Slope","Max","AUC","Stiffness", "Max", "Energy")
phen=c("BV.TV", "BV/TV")
pheno="Tb.Th3D"

p=df[[pheno]]
n=nrow(df)
sm=sum((p-mean(p))^2)
se=sqrt(sm/(n-1))

qqnorm((p-mean(p))/se)
qqnorm(log10(df[[pheno]]))

qqnorm(((df[[pheno]])-mean(df[[pheno]]))/sqrt((sd(df[[pheno]])))/)

qqline(log10(df[[pheno]]))
names(df)
names(biom2)
colnames(biom2)[6]="CC.Line"
df2=df
df$Mean2[108]=NA
df$Mean2=as.vector(df$Mean2)
setwd("c:/users/Roy/google drive/R scripts misc")
setwd("R scripts misc")
source("bars.R")
library("spdep")
library("Rcpp")
library("plotrix")
library("agricolae")
library("rafalib")
library("RColorBrewer")
pheno=c("BV.TV", 1)
p=df[[pheno]]
library(ggplot2)
unique(ds$group$M)



g1="IL-1513"
g2=c("IL-2126","IL-2391")
source("multiplot.R")


ph=c("BV.TV", "Tb.NPl", "Tb.Th","Conn.D", "BV/TV", "Tb.N", "Tb.Th","Conn.D")
ph=c("BMC", "Areal", "Length", "Cort.Th", "Mean2", "MaterProp", "Dia.Dia", "Polar")
ph=c("BVTV", "Tb.NPl", "Conn.D", "Tb.Th3D", "BV/TV", "Tb.N", "Conn.D", "Tb.Th")
ph=c("BV.TVP","Conn.DP","SMIP","BV.TV","Conn.D","SMIP","Tb.NP","Tb.ThP","SpP", "Mean1P", "BV.TVD", "Conn.DD", "SMID", "Tb.ND", "Tb.ThD", "SpD")
ph=c("BV.TVP","Conn.DP","SMIP","BV.TV","Conn.D","Tb.NP","Tb.ThP","SpP", "Mean1P", "BV.TVD", "Conn.DD", "SMID", "Tb.ND", "Tb.ThD", "SpD")

ph=c("vBMD","Porosity","BMC","Dia.Dia","Med.Dia","Cort.Th","Areal","BV.TV","Conn.D","SMI","Tb.N3D","Tb.NPl","Tb.Th", "Sp", "Mean1",
     "vBMD","Porosity","BMC","Dia.Dia","Med.Dia","Cort.Th","Areal","BV/TV","Conn.D","SMI","Tb.N3D","Tb.NPl","Tb.Th", "Sp", "Mean1")
ph=c("vBMD","Porosity","BMC")
ph=c("vBMD","Porosity")
ph=c("AUC","Slope","Max")
ph=phen
ph=c("BV.TV","Tb.NPl","Tb.Th","Conn.D","SMI","Sp", "BV/TV", "Tb.N", "Tb.Th", "Conn.D", "SMI", "Tb.Sp")
ph=c(phen, "BV/TV", "Conn.D", "Tb.N", "Tb.Th")
ph=c("vBMD","Cort.Th", "vBMD", "Ct.Th")
ph=c("Conn.D", "BV.TV")

phm=matrix(ph,2, byrow=T)
phm=matrix(ph,1, byrow=T)
pheno="vBMD"
pdf("PhenDistF.pdf")
p=list()

dfm=split(df, df$Sex)$M
dff=split(df, df$Sex)$F
par(mfrow=c(2,3))


df$BVTV=df$BVTV*100
df[["Mean2"]][108]=NA
df[["Cort.Th"]]=as.numeric(df[["Cort.Th"]])
i=7
i=1
df2=df[-which(df$CC.Line=="OR3393" | df$CC.Line=="OR15155"),]
df2=dfm
names(df)
p
df2=dff

df=dff
df2=biom2
df3=biom2
df3=df2

head(df2)
i=2

phm[1,]
df2=df
df3=df2
df2$Tb.NPl[-0]
p=list()
i=1
phm="vBMD"
df3=df
df2=df

df[[phm]][which(df[[phm]]=="#N/A")] = NA
for(i in 1:length(phm[1,])){

  #dfc=df[["Conn.D"]]
  #which(dfc=="#N/A")
  df[[phm[1,i]]][which(df[[phm[1,i]]]=="#N/A")] = NA
  if (sum(is.na(df2[[phm[1,i]]]))>0)
    df3=df2[-which(is.na(df2[[phm[1,i]]])),]
  fit=lm(as.numeric(df3[[phm[1,i]]])~CC.Line, data=df3)
  ds=duncan.test(fit, trt="CC.Line", console=F, alpha=0.001)
  dgrp=as.numeric(unique(ds$group$M))

  col=colorRampPalette(brewer.pal(10,"PuOr"))(length(dgrp))
  #source("bars.R"); 
  cat(phm[1,i], length(col), "\n",sep=" ")
  p[[phm[1,i]]]=descBars(col=col, grp=ds, 
                      df=df3, pheno=phm[1,i],ylb=phm[2,i],
                      desc=T, fill="#6699CC")+
                      theme(axis.text=element_text(size=15),
                      axis.title=element_text(size=15))
}

data.frame(ds$groups[,1])

p
ds

lapply()

df[which(df$CC.Line=="IL-519"),"vBMD"]


b = aggregate( df2$Cort.Th, by=list(CC.Line=df2$CC.Line), FUN="mean", na.rm=T)

se=sapply(df$CC.Line, function(x) std.error(df[["BV.TV"]]))

spdf2=split(df2, df2$CC.Line)
se=sapply(spdf2, function(x) std.error(x[["Cort.Th"]]))
data.frame(data.frame(table(df2$CC.Line))$Freq)


setwd("../YG/CC Proj/rhbdf")
df2=read.delim("BVTV36675S.txt")


b=data.frame(b)
b$se=se

b
df[[pheno]]=as.numeric(df[[pheno]])
mbvtvm=sapply(spdf2, function(x) mean(x[[pheno]]))

names(df)

p[[1]]

length(unique((df$CC.Line)))
sum(df$Sex=="F")
dfm=df[which(df$Sex=="M"),]
dff=df[which(df$Sex=="F"),]

length(unique(dfm$CC.Line))
length(unique(dff$CC.Line))

length(na.omit(df$Weight))
cor(na.omit(df$BV.TV)[1:160], na.omit(df$Mean1)[1:160])

summary(lm(na.omit(df$BV.TVD)[1:166]~na.omit(df$BV.TV)[1:166]))$r.squared

length(na.omit(df$Weight))
unique(df$CC.Line)
p
length(na.omit(df$Weight))
source("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/oct15/bone.R")

df[which(df$CC.Line=="IL-1513"),"BV.TV"]

bone.h2(file="BVTV36675.txt", col=phm[1,])

df[which(df$CC.Line=="IL-1513" & df$Sex=="M"),"BV.TV"]
p
sort(df$Tb.Th)
p
1104
p
plot(df2$BV.TV)
p[1]
p[2:length(p)]
ph[7]
sum(is.na(df[["Conn.DP"]]))
sum(is.na(df[["BV.TVP"]]))

names(df)
for (i in 1:4)
p[[i]]=p[[i]]+theme(text=element_text(size=2),axis.text.y=element_text(size=12), axis.text.x=element_text(size=3), axis.title.y=element_text(size=18))
p[[i]]=p[[i]]+geom_errorbar(size=5)

p[16]
source("c:/users/Roy/google drive/R scripts misc/multiplot.R")
setwd("c:/users/Roy/google drive/YG/thesis/")
getwd()
dev.off()
pdf("vertn.pdf")
 multiplot(p[[1]], NULL, NULL, p[[2]],NULL, NULL,cols=2)#,p[[5]],p[[6]], cols=2)
 multiplot(p[[1]], p[[2]], p[[3]], p[[4]], cols=2)
dev.off()            
,p[[7]],p[[8]],
            cols=2)
pdf("Figure 1.pdf")
  multiplot(p[[1]],p[[2]],p[[3]],p[[4]],cols=2)
p

dev.off()
names(biom2)

par(mfrow=c(2,2))
phm
i=1
for(i in 1:length(ph)){
 qqnorm(df[[phm[1,i]]],main=phm[2,i], col="red") 
 qqline(df[[phm[1,i]]])
}


i[[2]]
#+theme(legend.position="none")




r=order(p)
po=qnorm(r/(length(r) + 1))
sm=sum((p-mean(p))^2)
se=sqrt(sm/(n-1))
nr=(p-mean(p))/se

qqnorm(df[[pheno]])

names(df)
pdf("Trab_1192.pdf")
library(rafalib)
mypar2(1,2)
grp=ds$groups
sq[-c(g1,g2)]=2
df$grpBv=sq
agg=aggregate( df$grpBv, by=list(CC.Line=df$CC.Line), FUN="median" )
p + geom_segment(aes(fill=agg$x, x=0, y=mean(df[[pheno]]), xend=length(unique(df$CC.Line)), yend=mean(df[[pheno]])), lwd=1, lty=3)#, color="dark gray")+
  geom_segment(aes(x=0, y=mean(df[[pheno]])+sd(df[[pheno]]), xend=length(unique(df$CC.Line)), yend=mean(df[[pheno]])+sd(df[[pheno]])), lwd=0.5, lty=5)+#, color="dark gray")+
  geom_segment(aes(x=0, y=mean(df[[pheno]])-sd(df[[pheno]]), xend=length(unique(df$CC.Line)), yend=(mean(df[[pheno]])-sd(df[[pheno]]))), lwd=0.5, lty=5)+#, color="dark gray")+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))
dev.off()
length(unique(df$CC.Line))
p1=p+geom_segment(aes(x=2, 
                      y=0.25, xend=2, 
                      yend=0.35), size=1) +
  p1+geom_segment(aes(x=2, 
                   y=0.35, xend=28, 
                   yend=0.35), size=1) +

  p2=p1+geom_segment(aes(x=2, 
                      y=0.35, xend=28, 
                      yend=0.35), size=1) +
  p2+geom_segment(aes(x=28, 
                      y=0.35, xend=28, 
                      yend=0.15), size=1) 
  
b = aggregate( df$CC.Line, by=list(CC.Line=df$CC.Line), FUN="length")
choose(31,3)


#p+geom_vline(yintercept=5, color="green")

p1=p+geom_segment(aes(x=2.5, y=0.003, xend=8+.5, yend=0.003), color="blue", size=2)
p2=p1+geom_segment(aes(x=2.5, y=0.002, xend=25+.5, yend=0.002), color="dark blue", size=2)
p3=p2+geom_segment(aes(x=25.5, y=0.001, xend=31+.5, yend=0.001), color="dark red", size=2)
p4=p3+geom_line(data = arc.df, aes(Group+5.5, Value+0.07), lty = 2, lwd=1.0)

p5=p4+geom_line(data = arc.df2, aes(Group+2, Value+0.075), lty = 2, lwd=.7,color="red")
p6=p5+geom_segment(aes(x=2,xend=2,y=0.065,yend=0.07), lwd=.7, lty=2)
p7=p6+geom_segment(aes(x=9,xend=9,y=0.055,yend=0.07), lwd=.7, lty=2)
p8=p7+geom_segment(aes(x=1,xend=1,y=0.06,yend=0.0725), lwd=.7, lty=2)
p8+geom_segment(aes(x=1, y=0.0732, xend=3, yend=0.0732), lwd=.7, lty=2)
?geom_vline

t <- seq(0, 180, by = 1) / 180
t2<- seq(0, 1, by = .001) * pi / 180
r2 <- 10
x2 <- r2 * t
y2 <- 0.1*0.1 * t
length(y2)
y2=rep(0, 181)
seq
y2[20:162] <- y2[20] # Flattens the arc
length(t2)
arc.df2 <- data.frame(Group = x2, Value = y2)

r <- 3.5
x <- r * cos(t)
y <- 0.1*0.1 * sin(t)
y[20:162] <- y[20] # Flattens the arc

arc.df <- data.frame(Group = x, Value = y)



p3 + geom_path(aes(x=2.5, xend=25.5,y=0.04, yend=0.04))
names(df)
spl=split(df$Tb.Th3D, df$CC.Line)
t.test(spl[["IL-2680"]],rep(0.0484,5))

aov(Tb.NPl~CC.Line, data=df)

smfit=summary(lm(Conn.D~CC.Line, data=df))

fit=lm(df[[pheno]]~CC.Line, data=df)
summary(fit)
anova(fit)
AV=aov(fit)
library(agricolae)
ds=duncan.test(fit, trt="CC.Line", console=F, alpha=0.01)


grp=unique(ds$groups$M)
table(ds$groups$M,ds$groups$trt)
?duncan.test
grp1=sapply(grp, function(x) ds$groups$trt[which(ds$groups$M==x)])
grp1[[2]]

sp=split(df$Tb.Th3D,df$CC.Line)
x=sp[["IL-1061"]]
t.test(x,rep(mean(df$Tb.Th3D),3))
names(table(df$CC.Line))
geom(h=mean(df$Tb.Th3D))

(smfit)[4]
nm=names(which(smfit[4]))<0.001
pv=smfit[4][[1]][,4]
nm1=gsub("CC.Line", "",names(which(pv<=0.05 & pv>0.01)))[-1]
nm2=gsub("CC.Line", "",names(which(pv<=0.01)))[-1]
nm=gsub("CC.Line", "",nm)

se[nm]
lb.df=data.frame(l=nm, m=se[nm])


qqnorm(log10(df$Sp.dist))
names(df)
plot(df$weight, pch=16)
library(rafalib)
mypar2(2,1)
par(mfrow=c(3,2))


str(av)

pheno="Tb.N.prox"
st=(df[[pheno]]-mean(df[[pheno]]))/sd(df[[pheno]])
av=aov(st~CC.Line, data=df)
summary(lm(st~CC.Line, data=df))

label.df=data.frame(CC.Line=c("IL-1513","IL-2438"))

av
?aov
length(names(table(df$CC.Line)))
names(df)
(2*(1.96+1.03)^2/.7)*3.6
median(table(df$CC.Line))

2*(0.8+0.05)^2/105
2.4*1.5


setwd("R scripts misc")
scan=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scans/scanTh.txt")
dir()
ls()
(scan.data$scans$BVTV$genomewide.distribution)
scans=scan.data$scans$Conn.D
q = c( 0.5, 0.9, 0.95 )
quantiles = scans$genomewide.distribution[ floor( length(scans$genomewide.distribution)*q) ]
q
names(scans$result)
w+60:length(scans$result$BVTVP$logP)

mx2=which.max(scans$result$logP[w+2:length(scans$result$logP)])
w2=w+80
install.packages("qtl")
library("qtl")
data(hyper)
str(hyper)
(hyper$geno[[1]])
w=which.max(scan$logP)

scan$logP[w:(w+30)]
w2=w+31
w=1408
names(scans$result)
pos=(scan$from.bp[(w2-5):(w2+5)])
pos=as.numeric(as.character(pos))
lrt=(scan$logP[(w2-5):(w2+5)])
lrt=as.numeric(as.character(lrt))

qtl.ci<- function (pos = NULL, lrt = NULL, alpha = 0.05) {
  # pos: the vector of scanned positions
  # lrt: -2*log-likelihood ratio test statistic
  # alpha: the significance level
  tol<- alpha/100
  down<- 0.01 # drop extent each loop to search threshold
  np<- length(pos)
  # Step 1: Calculate RFR at each scanned position
  posq<- which.max(lrt)
  lrtq<- max(lrt)
  f<- exp((lrt - lrtq)/2)
  # Step 2: Calculate the area Ai for position i to i+1
  Area<- (pos[2:np] - pos[1:(np-1)]) * (f[1:(np-1)] + f[2:np])/2
  # Step 3: Scale Ai
  SArea<- sum(Area)
  Area<- Area/SArea
  # Step 4: calculate CI
  thr<-lrtq
  loop<-TRUE
  rp<- rev(pos)
  rf<- rev(f)
  while(loop)
  {
    thr<-thr-down
    fd<-exp((thr-lrtq)/2)
    tArea<-lArea<-rArea<-dl<-dr<-0
    u<-(lrt<=thr)
    v<-rev(u)
    u<- sum(cumprod(u))
    v<- sum(cumprod(v))
    if (u>0)
    {
      dl<- (pos[u+1]-pos[u])*(fd - f[u])/(f[u+1]-f[u])
      lArea<- 0.5*dl*(f[u] + fd)/SArea
    }
    if (u>1) lArea<- lArea + sum(Area[1:(u-1)])
    if (v>0)
    {
      dr<- (rp[v]-rp[v+1])*(fd - rf[v])/(rf[v+1]-rf[v])
      rArea<- 0.5*dr*(rf[v] + fd)/SArea
    }
    if (v>1) rArea<- rArea + sum(rev(Area)[1:(v-1)])
    tArea<- lArea + rArea
    if ((tArea<=alpha)|(abs(tArea-alpha)<=tol)) loop<-FALSE
  }
  if (u == 0) u<-1
  if (v == 0) v<-1
  lower<- pos[u] + dl
  upper<- rev(pos)[v] - dr
  return(c(lower,upper))
}
#End
getwd()
source("R scripts misc/qtlcsi.R")


ci=qtl.ci(pos=as.numeric(pos), lrt=as.numeric(lrt), alpha = c(0.5))
ci/1e6

min(pos)
(ci[2] - ci [1])/1e6
which(pos==1031)
pos
w=which.max(scans$result$logP)
scans$result$chromosome[w+100]
gt=(scans$result$genomewide.permute.pval)
res=p.adjust(scans$result$logP)
range(res)
s=seq(w2-50,w2+50)
range(pos)


set.seed(123)
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))

pc=p.adjust(p)

s

sum(p2<0.01)

list(p2)
write.table(p2,file="p2.txt",sep="\t",quote=FALSE,row=FALSE)

p2=10^(-lrt)
pa=p.adjust(p2, "BH")

sum(pa<0.000005)

bf=0.05/length(p2)
sum(p2<=bf)
library("fdrtool")

hc=hc.score(p2)
hc.thresh(p2)



s=rep(1,1000)
s[950:1000]=0.0000001

padj2=p.adjust(s, "BH")
sum(padj2<0.00001)

s1=sum(lrt>4)
pa=p.adjust(10^(-lrt), "BH")
s2=sum(pa<0.05)

s1/s2

sum(p<0.05)
padj=p.adjust(p, "BH")
sum(padj<0.05)
plot
matplot(p2, pa, type="l", asp=1)

p <- 2*pnorm(sort(-abs(x)))
## or all of them at once (dropping the "fdr" alias):
p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
p.adj    <- sapply(p.adjust.M, function(meth) p.adjust(p, meth))
matplot(p2, p.adj, ylab="p.adjust(p, meth)", type = "l", asp = 1, lty = 1:6,
        main = "P-value adjustments")
legend(0.7, 0.6, p.adjust.M, col = 1:6, lty = 1:6)


pos=as.numeric(as.character(scans$result$from.bp[s]))
lrt=as.numeric(as.character(scans$result$logP[s]))
lrt2=-2*log(10^(-lrt))
which.max(lrt2)
scans$result$from.bp[s[21]]
hist(f)
hist(log10(lrt))
lrt3=(log10(lrt))
rn=range(pos)
plot(f, type="l")
plot(Area)
length(pos)
rn[2]-rn[1]
ci[2] - ci [1]
length(scans$result$from.bp[w-100:w+100])
s=seq(w-100,w+100)
w+100
pos=as.numeric(as.character(scans$result$from.bp[s]))
range(pos)
w=as.numeric(w)
class(1)

z=c(1,2,3,4,5)
z[1:3]





#########heritability may 1st 2015
setwd("../YG/CC Proj")

getwd()
setwd("redHawk/scripts")
source("heritability.R")

df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scan_1192.txt")

file="http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scan_TbTh601.txt"
file="http://mus.well.ox.ac.uk/preCC/BONE/pheno.txt"
df=read.delim(file)
file="http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/march15/imputed15.txt"
file="http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/prelim103.txt"
length(names(table(df$CC.Line)))


col=c(17:43)
trab.h2(file=file, txt.file="OldPhenoH2.txt", col=col)

names(df)

sink()
dev.off()


cond=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/Conn.D.scan.txt")
names(cond)
gw=list()
gw$SNP=bv$from.marker
gw$CHR=bv$chromosome
gw$BP=bv$from.bp
gw$P=10^(-bv$logP)

gw$CHR=(gsub("chr","",gw$CHR))

gw$CHR[which(gw$CHR=="X")]="20"
gw$CHR=as.numeric(gw$CHR)
manhattan(gw, cex=.5)

plot(cond$logP, type="l")
unique(gwasResults$CHR)
nrow(data.frame(gw))
install.packages("qqman")
library(qqman)
qq(gw$P)
qqnorm(gw$P)
bv=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/BVTV.scan.txt")
ls()
remove(cd)
c12=cond[which(cond$chromosome=="chr12"),]

10^(-bv$logP[order(bv$logP, decreasing=T)[1:5]])

smbv=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim.txt"); plot(smbv[,1][which(smbv[,1]>6)], type="b", pch=19)
sm=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_ConnD1.txt"); plot(sm[,1][which(sm[,1]>5.5)], type="b", pch=19)
sm2=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_TbTh.txt"); plot(sm2[,1][which(sm2[,1]>6)], type="b", pch=19)
sm3=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_ConnD12.txt"); plot(sm3[,1][which(sm3[,1]>5.9)], type="b", pch=19)
sm4=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_ConnD.14.txt"); plot(sm4[,1][which(sm4[,1]>5.9)], type="b", pch=19)
sm5=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_TbNPl.19.txt"); plot(sm5[,1][which(sm5[,1]>5.5)], type="b", pch=19)
sm6=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/MCSim_TbN.txt"); plot(sm6[,1][which(sm6[,1]>6)], type="b", pch=19)
nrow(smbv)
th=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/Tb.Th3D.scan.txt")
plot(th$logP, type="l")
sm$V1
sam=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/linesRemoved1638.txt")
sam
0
cb=(0.7*bv$logP+0.3*cond$logP)
par(mfrow=c(1,1))
plot(cond$logP, type="l")



th=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/Tb.Th3D.scan.txt")
lr=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/linesRemoved_TbTh601.txt")
c7=th[which(th$chromosome=="chr7"),]
max(c7$logP)

setwd("YG/CC Proj/redHawk/")
getwd()
setwd("../scripts")

load("dbMatTbTh601.RData")
dim(d.1)
load("phen.dataTbTh601.RData")
load("mat185_36TbTh.RData")
load("matRand.RData")
ls()
mt=mat
mr
####Correlation Matrix
      library("sna")
      #db=condensed.db
      #df=read.delim("imputed15.txt")
      #d.1 <- db$additive$chr$chr8[156][[1]]$mat
      
      nm=db$subjects[which(db$subjects %in% names(table(df$CC.Line)))]
      nm=(phen.data$CC.Line)
      
      mtx=matrix(nrow=length(nm), ncol=length(nm))
      mt=d.1[which(db$subjects %in% names(table(df$CC.Line))),]
      mt=d.1
      mt=mr
      diag(mtx)=1
      for(i in 1:length(nm)) {
        j=i+1
        while(j<=length(nm)) {
          mtx[j,i]=cor(mt[i,],mt[j,])
          j=j+1
        }
      }

head(mtx)
mtxs=symmetrize(mtx,"lower")
mtxs2=symmetrize(mtx,"lower")
setwd("../YG/CC Proj/redHawk/scripts")
load("mtx.RData")
load("mtxs.RData")

int=colnames(mtxs)

df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/march15/imputed15.txt")
df=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/prelim103.txt")
agg.phen.data = aggregate( df[["BV.TV.met"]], by=list(CC.Line=df$CC.Line), FUN="median" )

sms=split(df, df$Sex)
sms$M
par(mfrow=c(1,1))
layout(matrix(c(1,2,3,3,4,5), 2, 3, byrow = TRUE))

?layout


pheno="BVTV"
boxplot(sms$M[[pheno]], sms$F[[pheno]], notch=T, range=1, col=c("#99CCFF", "#CC66FF"), main=pheno, names=c("Male","Female"))

#which(agg.phen.data$CC.Line %in% int)
phen.data = agg.phen.data[agg.phen.data$CC.Line %in% int,]

mtxs
plot(mtxs)


pv=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/Tb.N3D.scan.txt")
pvals=10^(-pv$logP)
c7=pv[which(pv$chromosome=="chr1"),]
pvals=10^(-c7$logP)
  
qq(pvals, cex=0.5)
mrg=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/BVTV.merge/BVTV.chr4:16200000-17900000.merge.txt", sep=" ")

names(mrg)      
pv=mrg$logP.merge
pvals=10^(-pv)
      observed <- sort(pvals); lobs <- -(log10(observed)); expected <- c(1:length(observed)); lexp <- -(log10(expected / (length(expected)+1)))
      plot(lexp, lobs, pch=23, cex=.4, bg="black") 
      plot(lobs, lexp, pch=23, cex=.4, bg="black") 
      plot(lobs, cex=0.5) 
      points(lexp, col="red", lwd=0.5, cex=0.5) 
      abline(0,1, col="red")
      
      pdf("qqplot.pdf", width=6, height=6)
      plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
      
      dev.off()




library("qqman")


qq(10^(-pv$logP), type="p", pch=1)


library(gplots)
heatmap(mtxs)
qq
pvector=pvals
o = -log10(sort(pvector, decreasing = FALSE))
n=(length(pvector))
e = -log10(ppoints(n, a = ifelse(n <= 10, 3/8, 1/2)))
e=-log10(sort(seq(1:n)/n, decreasing=T))
abline(0,1, col="red")
qqnorm(e)

plot(e)
points(o, col="red")
plot(x=e, y=o)#,pch = 20, xlim = c(0, max(e)), ylim = c(0,
                                             max(o)), xlab = expression(Expected ~ ~-log[10](italic(p))),
     ylab = expression(Observed ~ ~-log[10](italic(p))))

ppoints
library("corrplot")
corrplot(mtxs, cl.pos="b", type="lower", col=hmcol,tl.cex=0.5, method="square", tl.col="black")
?corrplot
#3.week2.2.R
library("RColorBrewer")
library(rafalib)
hmcol=colorRampPalette(brewer.pal(9,"GnBu"))(100)
?RColorBrewer
mypar2(1,1)
?brewer.pal
setwd("../../../../R scripts misc")
pdf("kinshipPlotTbThChr7.pdf")

image(mtxs2)
heatmap(mtxs2)#,col=hmcol)
dev.off()
heatmap.2(e[idx,],labCol=tissue,trace="none", ColSideColors=cols,col=hmcol)

library(lattice)
levelplot(mtxs2)#, col.regions= hmcol)
cols2 <- colorRampPalette(c("blue","white","red"))(256)
image(mtxs, col=hmcol)

L=t(mtxs2)
?chol
chol(mtxs)
library("PropCIs")
mtxs2==mtxs
mtxs2=sweep(mtxs2, 2, colSums(mtxs2), FUN="/")
mtxs2=scale(mtxs2, center=FALSE, scale=colSums(mtxs2))
chol(mtxs2)
mtxs2=t(apply(mtxs, 1, function(x)(x-min(x))/(max(x)-min(x))))

chol(mtxs2)
library("regress")



r=regress(phen.data$x~1,~mtxs)
fit0=lm(phen.data$x~1)
fit0=lm(rn~)
fit00=lm(phen.data$x~sex)
fit000=lm(phen.data$x~sex+mtxs3)
fit4=lm(phen.data$x~d.1)#mtxs[,10]+mtxs[,2]+mtxs[,3]+mtxs[,4]+
          mtxs[,5]+mtxs[,6]+mtxs[,7]+mtxs[,8]+
          mtxs[,8]+mtxs[,9]+mtxs[,9]+mtxs[,10])
class(nm)
mtxs[,10]
#fit4=lm(rn~d.1)
summary(fit4)

rn=rnorm(39)

a=anova(fit0,fit00, fit000)
a=anova(fit0,fit4, test="Chisq")
?anova


h = (a[2,2]-a[3,2])/(a[2,2])

sigma=r$sigma
sigma.cov=r$sigma.cov
install.packages("gap")
library(gap)
h2=h2G(sigma,sigma.cov)
install.packages("hglm")
library(hglm)
phen.data$CC.Line
sx=table(df$CC.Line,df$Sex)
sx[which(sx>1)]=1
sx[which(sx<1)]=0

model.matrix(~sex)
dim(sx)

sex=sx[which(rownames(sx) %in% phen.data$CC.Line),]

library(hglm)
library(ggplot2)

unloadNamespace("ggplot2")
unloadNamespace("Hmisc")
unloadNamespace("MASS")
unloadNamespace("hglm.data")
unloadNamespace("biovizBase")
unloadNamespace("Gviz")
unloadNamespace("interactiveDisplay")
unloadNamespace("AnnotationHub")

sxd=rep(1,14)
sxd2=rep(0,14)
sxd=c(sxd,sxd2)

df$Sex

?hglm
m1=hglm(phen.data$x, X=model.matrix(~sxd), Z=mtxs2)
m12
summary(m1)
m1$varRanef/(m1$varRanef+m1$varFix)
x=phen.data$x
sm=sum((x-mean(x))^2)
se=sqrt(sm/(length(x)-1))
x=((x-mean(x))/se)
nmt=sapply(1:39, function(i) {
  sm=sum((mtxs[,i]-mean(x))^2)
  se=sqrt(sm/(length(mtxs[,i])-1))
  ((mtxs[,i]-mean(mtxs[,i]))/se)
})  


diag(nmt)=1
sex[,1][which(sex[,1]>0)]=1
sex[,2][which(sex[,2]>0)]=1

m1=hglm(y=x, X=model.matrix(~sex), Z=mtxs)

?hglm
m1$varRanef/(m1$varRanef+m1$varFix)

mtxs3=mtxs
mtxs3[which(mtxs3>0.35 & mtxs3<0.65)]=0.5
mtxs3[which(mtxs3>.15 & mtxs3<=0.35)]=0.25
mtxs3[which(mtxs3>0 & mtxs3<=.15)]=0.125
mtxs3[which(mtxs3>0.35 & mtxs3<0.4)]=0.25
mtxs3[which(mtxs3>=0.4 & mtxs3<0.77)]=0.5
mtxs3[which(mtxs3>=.77 & mtxs3<1)]=0.5
mtxs3[which(mtxs3<0)]=0
diag(mtxs3)=1





##########broad regional
setwd("YG/CC Proj/redHawk/scripts")

tb=read.delim("http://www.broadinstitute.org/files/shared/diabetes/scandinavs/DGI_chr9_hit.txt")
tb2=tb[,-5]
names(tb)
tb2$RSQR=0
tb2$RSQR[which(tb2$TYPE=="typed")]=tb$RSQR[which(tb2$TYPE=="typed")]
tb2$RSQR[50:70]=tb$RSQR[50:70]
setwd("../scripts")
source("fancyAssocplot.R")
make.fancy.locus.plot("rs10811661", "CDKN2A/CDKN2B region", "9", tb2, 9, 5.4e-8)

kg=read.delim("known_genes_build35_050307_chr1.txt", sep=" ")
names(kg)
tb2["rs10811661",]

kg[1,]
head(tb2)


mrg4=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/mf/7jan15/trab/BVTV.merge/BVTV.chr4:16200000-17900000.merge.txt", sep=" ")
mrg2=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/BVTV.merge/BVTV.chr2:1.22e+08-130700000.merge.txt", sep=" ")
names(mrg4)
bv=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/BVTV.scan.txt")
c4=bv[which(bv$chromosome=="chr2"),]
c4p=c4[which(c4$from.bp>=130400000 & c4$to.bp<=130700000),]
1.22e+08-130700000
tbx=list()

tbx$SNP=c4p$from.marker
tbx$POS=c4p$from.bp
tbx$PVAL=10^(-c4p$logP)
tbx$TYPE=rep("typed")
tbx$RSQR=0
tbx=data.frame(tbx)
head(tbx)
tb2
mrgp=mrg2[which(mrg2$logP.merge>3),]
mrgp=mrgp[which(mrgp$bp>130550000 & mrgp$bp<130590000),]
names(mrgp)

make.fancy.locus.plot("JAX00546093", "Test region", "4", tbx, 8, pmin)
pmin=min(tbx$PVAL)
snp="JAX00546093"

rec=read.delim("known_genes_build35_050307_chr4.txt", sep=" ")
head(rec)
gm=read.delim("genetic_map_chr4.txt", sep=" ")
max(gm2[,2])
gmx=list()
gm2=gm
head(tb)
gm2[,3]=0
names(gm)
gmx$SNP=factor(seq(500:(500+nrow(mrgp)))[-1])
gmx$POS=mrgp$bp
gmx$PVAL=10^(-mrgp$logP.merge)
gmx$TYPE=rep("imputed", nrow(mrgp))
gmx$RSQR=rep(0, nrow(mrgp))
length(gmx$SNP)
gmx
names(mrgp)
gmx$COMBINED_rate.cM.Mb.=mrgp$logP.merge
gmx$Genetic_Map.cM.=0
levels(tbx[,1])=c(levels(tbx[,1]), gmx$SNP)
gmx=data.frame(gmx)
dfx=rbind(tbx, gmx)
head(dfx)
unique(gmx[,2])
locus=tbx
locus[snp,]
locus

snp=c4$from.marker[which.max(c4$logP)]
dfx$RSQR[which(dfx$TYPE=="typed")]=0.5
dfx$RSQR[which(dfx$TYPE=="imputed")]=0.2
make.fancy.locus.plot(snp, "Test region", "2", gmx, 8, pmin)

kg=read.delim("knownGene.txt", header=F)



###cisim
load("YG/CC Proj/redHawk/scripts/cisimdf2.RData")

(data.frame(df2))

f1=as.numeric(sapply((df2[1,]), function(x) x))
f2=as.numeric(sapply((df2[4,]), function(x) x))
f3=as.numeric(sapply((df2[5,]), function(x) x))
plot(ecdf(f2), lty=5, add=T)

hist((f3/1e+6))
boxplot(f1~f2, range=0.5)
hist(sapply(df2[4,], function(x) x))
boxplot(tdf2[,1]~tdf2[4,], range=1)

phen="Tb.Th.met"
h=hist(df[[phen]])
h$density = h$counts/sum(h$counts)
dens=density(df[[phen]])
plot(h, freq=F, ylim=c(0,1))
lines(ecdf(df[[phen]]))

points(dens$x, dens$y/100, col="blue", type="l",lwd=0.5)
str(dens)
?density

setwd("../YG/CC Proj/redHawk/avp")

load("../scripts/cisimTh.RData")
dfth

f1th=as.numeric(sapply((dfth[1,]), function(x) x))
f2th=as.numeric(sapply((dfth[2,]), function(x) x))
f3th=as.numeric(sapply((dfth[3,]), function(x) x))
f4th=as.numeric(sapply((dfth[4,]), function(x) x))
f5th=as.numeric(sapply((dfth[5,]), function(x) x))

h=hist(f5th)
h$density = h$counts/sum(h$counts)
h$breaks
hist(f4th, breaks=30, xlim=c(-10,10), col="gray",freq=F)
plot(h, freq=F, col="gray")
boxplot(f3th, f2th)
plot(ecdf(f1th))#, xlim=c(0,12))
lines(ecdf(f2th), col="red")
hist(dfth[[1,14]])
plot(dfth[[4,1]],dfth[[1,1]])
rownames(dfth)

##FDR

#bv=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scans/scanBV.txt")
#th=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scans/scanTh.txt")
#tb=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scans/scanTn.txt")
#con=read.delim("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/scans/scanCon.txt")


?p.adjust
arr=list(bv, tb,th,con)
for (ph in arr) {
    df=ph
    #print(df)
    pv=10^(-df$logP)
    pa=p.adjust(pv, "BH")
    
    idx=which(df$logP>5)
    
    cat("FDR =",100*(1-length(which(pa[idx]<0.01))/length(idx)),"%","\n")
}

install.packages("gdata")
library("gdata")
read.xls
setwd("../../../batch sheets")
am=read.csv("all march 15.csv")
head(am)
am2=read.csv("21st batch suggested.csv")
head(am2)
t2=(table(am2[[2]],am2[[5]]))
colnames(t1)
t4=rbind(t1,t2)

rownames(t3)
t3=data.frame(t3)
t3[[1]]=rownames(t4)




length(am2[[2]])
am2=am2[1:34,1:7]

am=am[1:293,1:7]
colnames(am2)=c(1:7)
am3=rbind(am,am2)
am4=am3
for(i in 1:nrow(am4)){
  if (am4[i,7]=="B" && am4[i,2]<=400) am4[i,2]=(am4[i,2]+500)
}
am4[[2]]

am3[[2]]+500
t2=table(am3[[2]],am3[[5]])
t3[,1]=rn
rn2=rownames(t2)

length(rn2)


write.table(t3,file="jun15Summary.txt",sep="\t",quote=FALSE,row=FALSE)


source("http://mus.well.ox.ac.uk/preCC/BONE/LEVY/redHawk/avp/CIfigs.R")
setwd("C:/YG/CC Proj/rhbdf")
data.df=read.delim("cisim/simCI.11.536.Tb.NPl.txt", sep="")
calc.CIs(data.df, low=0.1,high=0.9)







#march 2016
setwd("c:/users/Roy/google drive/YG/CC Proj/rhbdf/")
rd=read.delim("rhbdftab4.txt")
phen="SMI"
lm1=lm(rd[[phen]]~1, data=rd)
lm2=lm(rd[[phen]]~strain+weight+age+batch, data=rd)
an=(anova(lm1,lm2))$Pr[2]



nm=c(names(rd)[c(5:12)],"BV/TV","Tb.N","Tb.Th","Conn.D","SMI","Tb.Sp","vBMD","Ct.Th")

names(rd)
pdf("Fig.7n2.pdf")
par(mfrow=c(3,3))
i=9
for(phen in nm[1:8]){
  wna=which(is.na(rd[[phen]]))
  if(sum(wna)>0) rdx=rd[-wna,]
  else rdx=rd
  boxplot(rdx[[phen]]~strain, data=rdx, main=nm[i], notch=F, col="darkblue", frame.plot=F)
  lm1=lm(rdx[[phen]]~1, data=rdx)
  lm2=lm(rdx[[phen]]~strain+weight+age+batch+season+year, data=rdx)
  #an=(t.test(rdx[[phen]][which(rdx$strain=="ko")], rdx[[phen]][which(rdx$strain=="wt")]))$p.value
  #may292016
  pvn=t.test(rd[(rd$strain=="KO"),phen],rd[(rd$strain=="WT"),phen])$p.value
  an=(anova(lm1,lm2))$Pr[2]
  subtitle=paste("PV =",round(an,3), "(",round(pvn,3),")" )  
  mtext(subtitle, cex=0.8)
  stripchart(rdx[[phen]]~strain, data=rdx, add=T, vertical=T, method="jitter", col="blue", pch=19, lwd=3)
  i=i+1
}
dev.off()

sum(rd$strain=="wt")
