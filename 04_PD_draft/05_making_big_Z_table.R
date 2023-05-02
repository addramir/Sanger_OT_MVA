library("data.table")

setwd("~/projects/MVA_output/02_PD_gwas/full_gwas/")
#load("~/projects/MVA_output/02_PD_gwas/v1_hm3_4x4_matrix.RData")
lst=c("nallsEtAl2019_excluding23andMe.txt","NEALE2_20107_11.txt","NEALE2_20110_11.txt",
"FINNGEN_R6_G6_PARKINSON_INCLAVO.txt","FINNGEN_R6_G6_PARKINSON.txt")


i=1
fl=lst[i]
x=fread(fl,data.table=T)
x=x[,Z:=b/se]
x=x[,c("SNP","Z","N","eaf")]
colnames(x)=c("SNP",paste(c("Z","N","eaf"),fl,sep="_"))
out=x
out=out[SNP!=""]


i=2
for (i in 2:length(lst)){
	fl=lst[i]
	x=fread(fl,data.table=T)
	x=x[,Z:=b/se]
	x=x[,c("SNP","Z","N","eaf")]
	x=x[SNP!=""]
	colnames(x)=c("SNP",paste(c("Z","N","eaf"),fl,sep="_"))
	out=merge(out,x,by="SNP",all=T)
}

save(list=c("out","lst"),file="~/projects/MVA_output/02_PD_gwas/5_Traits_Z_matrix.RData")

i=1
j=1
n=length(lst)
phe=array(NA,c(n,n))
diag(phe)=1
maf_thr=0.05
for (i in 1:(n-1)){
	for (j in (i+1):n){
		ind=c((i-1)*3+2,(j-1)*3+2)
		Z=out[,ind,with = FALSE]
		Z=as.data.frame(Z)
		ind=c((i-1)*3+4,(j-1)*3+4)
		EAF=out[,ind,with = FALSE]
		EAF=as.data.frame(EAF)

		EAF[EAF<maf_thr | EAF>(1-maf_thr)]=NA
		Z[Z^2>4]=NA

		Z=cbind(Z,EAF)
		Z=na.omit(Z)

		phe[i,j]=phe[j,i]=cor(Z[,1],Z[,2])
	}
}

save(list=c("out","lst","phe"),file="~/projects/MVA_output/02_PD_gwas/5_Traits_Z_matrix.RData")

