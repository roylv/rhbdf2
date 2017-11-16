for (file in c("alBV.TVSx.txt",
"alvBMDSx.txt",
"alCort.ThSx.txt",
"alSpSx.txt",
"alSMISx.txt",
"alConn.DSx.txt",
"alTb.ThSx.txt",
"alTb.NPlSx.txt")) {
	df=read.delim(file)
	df=df[df$Sex=="M",]; df=df[(df[[3]] %in% c("A","C")),]	
	cat("\n",file,"\n")
	print(t.test(df[df[[3]]=="A",1], df[df[[3]]=="C",1]))


}

