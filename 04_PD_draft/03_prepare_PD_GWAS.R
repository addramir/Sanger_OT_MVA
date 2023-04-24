library("data.table")

setwd("~/projects/MVA_output/02_PD_gwas/gwas/")

ref=fread("NEALE2_20107_11.txt",data.table=F)

lexa1a2=function(y){
	out=apply(y,MAR=1,function(x){paste(sort(x),collapse=":")})
	out
}

out=lexa1a2(ref[,c("A1","A2")])
snpid_ref=paste(paste("chr",ref$chr,sep=""),ref$pos,out,sep=":")
ref=cbind(ref,snpid=snpid_ref)

####
fl=fread("../PD_gwas/nallsEtAl2019_excluding23andMe_allVariants.tab",data.table=F)
out=lexa1a2(fl[,c("A1","A2")])
snpid_fl=paste(fl$SNP,out,sep=":")

table(snpid_ref%in%snpid_fl)

ind=which(snpid_fl%in%snpid_ref)
fl=fl[ind,]
snpid_fl=snpid_fl[ind]
ind=match(snpid_fl,snpid_ref)
table(snpid_fl==snpid_ref[ind])
fl=cbind(fl,rsid=ref[ind,"SNP"])
fl=fl[!(fl$p==1),]
fl=fl[fl$b!=0,]

fls=fl[,c("rsid","A1","A2","b","se")]
fls=cbind(fls,N=fl$N_cases+fl$N_controls)
colnames(fls)=c("SNP","A1","A2","b","se","N")
fwrite(file="nallsEtAl2019_excluding23andMe.txt",x=fls,sep="\t")

####
fl=fread("../PD_gwas/IPDGC_AAO_GWAS_sumstats_april_2018.txt",data.table=F)
out=lexa1a2(fl[,c("Allele1","Allele2")])
snpid_fl=paste(fl$MarkerName,toupper(out),sep=":")


table(snpid_ref%in%snpid_fl)

ind=which(snpid_fl%in%snpid_ref)
fl=fl[ind,]
snpid_fl=snpid_fl[ind]
ind=match(snpid_fl,snpid_ref)
table(snpid_fl==snpid_ref[ind])
fl=cbind(fl,rsid=ref[ind,"SNP"])
fl=fl[fl[,"P-value"]!=1,]
fl=fl[fl[,"Effect"]!=0,]

fls=fl[,c("rsid","Allele1","Allele2","Effect","StdErr")]
#fls=cbind(fls,N=fl$N_cases+fl$N_controls)
colnames(fls)=c("SNP","A1","A2","b","se")
fls=cbind(fls,N=17996)
fls[,"A1"]=toupper(fls[,"A1"])
fls[,"A2"]=toupper(fls[,"A2"])
fwrite(file="IPDGC_AAO_GWAS_sumstats_april_2018.txt",x=fls,sep="\t")

####
fl=fread("../PD_gwas/PDsubtype_Binary_summary_stats.txt",data.table=F)
out=lexa1a2(fl[,c("Allele1","Allele2")])
snpid_fl=paste(fl$MarkerName,toupper(out),sep=":")

table(snpid_ref%in%snpid_fl)

ind=which(snpid_fl%in%snpid_ref)
fl=fl[ind,]
snpid_fl=snpid_fl[ind]
ind=match(snpid_fl,snpid_ref)
table(snpid_fl==snpid_ref[ind])
fl=cbind(fl,rsid=ref[ind,"SNP"])
fl=fl[fl[,"P-value"]!=1,]
fl=fl[fl[,"Effect"]!=0,]

fls=fl[,c("rsid","Allele1","Allele2","Effect","StdErr")]
#fls=cbind(fls,N=fl$N_cases+fl$N_controls)
colnames(fls)=c("SNP","A1","A2","b","se")
fls=cbind(fls,N=3212)
fls[,"A1"]=toupper(fls[,"A1"])
fls[,"A2"]=toupper(fls[,"A2"])
fwrite(file="PDsubtype_Binary_summary_stats.txt",x=fls,sep="\t")
