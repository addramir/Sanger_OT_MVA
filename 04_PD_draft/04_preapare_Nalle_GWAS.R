library("data.table")

setwd("~/projects/MVA_output/02_PD_gwas/full_gwas/")

ref=fread("NEALE2_20107_11.txt",data.table=F)

lexa1a2=function(y){
	y1=y[,1]
	y2=y[,2]
	ind=which(y1>y2)
	y3=y1
	y1[ind]=y2[ind]
	y2[ind]=y3[ind]
	out=paste(y1,y2,sep=":")
	#out=apply(y,MAR=1,function(x){paste(sort(x),collapse=":")})
	out
}

ref=ref[ref$SNP!="",]

out=lexa1a2(ref[,c("A1","A2")])
snpid_ref=paste(paste("chr",ref$chr,sep=""),ref$pos,out,sep=":")
ref=cbind(ref,snpid=snpid_ref)



####
fl=fread("../PD_gwas/nallsEtAl2019_excluding23andMe_allVariants.tab",data.table=F)
out=lexa1a2(fl[,c("A1","A2")])
snpid_fl=paste(fl$SNP,out,sep=":")

table(snpid_fl%in%snpid_ref)

ind=which(snpid_fl%in%snpid_ref)
fl=fl[ind,]
snpid_fl=snpid_fl[ind]


ind=match(snpid_fl,snpid_ref)
table(snpid_fl==snpid_ref[ind])
fl=cbind(fl,rsid=ref[ind,"SNP"])

ref=ref[ind,]

ind=which(fl$A1!=ref$A1)
fl$b[ind]=-fl$b[ind]
fl$freq[ind]=1-fl$freq[ind]

fl$A1=ref$A1
fl$A2=ref$A2


fl=fl[!(fl$p==1),]
fl=fl[fl$b!=0,]


fls=fl[,c("rsid","A1","A2","b","se","freq","p")]
fls=cbind(fls,N=fl$N_cases+fl$N_controls)
colnames(fls)=c("SNP","A1","A2","b","se","eaf","pval","N")
fwrite(file="nallsEtAl2019_excluding23andMe.txt",x=fls,sep="\t")




LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_hapmap2_SVD_eigen99_extraction"

gwas1=fread("~/projects/MVA_output/02_PD_gwas/gwas/nallsEtAl2019_excluding23andMe.txt",data.table=F)
gwas2=fls
library("HDL")
res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)