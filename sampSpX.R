setwd("/scratch200/roeilevi/rhbdf")

source("happyEx19.preCC.R")

condensed.db = load.condensed.database("../Mar2015.MEGA+MDA+MUGA/CONDENSED/")

source("~/public_html/rhbdf/happyMOD.R")

scan.phenotypes.lm2(smn=5, phn="SpX",chrMax="chrX", rand=1, condensed.db=condensed.db, data.file="~/BVTV36675.txt", phen="Sp", permute=150, covariates = c("Sex", "WeekOld","Batch" ), aggregate.data="residuals", use.weights=T, normalise=FALSE, zero.na=FALSE, cousins=TRUE, scan.file="scans_SpX.RData", mc.cores=40,n=1,mast="FemOctNoOut.txt", AddMast=4)

















