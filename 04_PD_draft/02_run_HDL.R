library("data.table")
library("HDL")

setwd("~/projects/MVA_output/02_PD_gwas/gwas/")
lst=list.files()
n=length(lst)

LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_SVD_eigen99_extraction"
h2=rep(NA,n)
h2_se=rep(NA,n)
h2_p=rep(NA,n)
intspt=rep(NA,n)
for (i in 1:n){
	gwas1=fread(lst[i],data.table=F)
	res.HDL <- HDL.h2(gwas1, LD.path,intercept.output=TRUE)
	h2[i]=res.HDL$h2
	h2_se[i]=res.HDL$h2.se
	h2_p[i]=res.HDL$P
	intspt[i]=res.HDL$intercept
}

LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_hapmap2_SVD_eigen99_extraction"
h2_hm2=rep(NA,n)
h2_hm2_se=rep(NA,n)
h2_hm2_p=rep(NA,n)
intspt_hm2=rep(NA,n)
for (i in 1:n){
	gwas1=fread(lst[i],data.table=F)
	res.HDL <- HDL.h2(gwas1, LD.path,intercept.output=TRUE)
	h2_hm2[i]=res.HDL$h2
	h2_hm2_se[i]=res.HDL$h2.se
	h2_hm2_p[i]=res.HDL$P
	intspt_hm2[i]=res.HDL$intercept
}

save(list=c("lst","h2","h2_se","h2_p","intspt"),file="~/projects/MVA_output/02_PD_gwas/v1_hm3_h2.RData")
save(list=c("lst","h2_hm2","h2_hm2_se","h2_hm2_p","intspt_hm2"),file="~/projects/MVA_output/02_PD_gwas/v1_hm2_h2.RData")

#############

ind=which(h2_p<=0.05/20)
lst=lst[ind]
n=length(lst)
LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_SVD_eigen99_extraction"
out_rg=array(NA,c(n,n))
out_pval=array(NA,c(n,n))
for (i in 2:(n-1)){
	for (j in (i+1):n){
		idd=paste(i,j,sep=" ")
		if (!(idd %in%c("2 4"))){
			gwas1=fread(lst[i],data.table=F)
			gwas2=fread(lst[j],data.table=F)
			res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)
			out_rg[i,j]=out_rg[j,i]=res.HDL$rg
			out_pval[i,j]=out_pval[j,i]=res.HDL$P
		}	
	}
}

save(list=c("lst","out_rg","out_pval","h2"),file="~/projects/MVA_output/02_PD_gwas/v1_hm3_4x4_matrix.RData")

######

#LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_hapmap2_SVD_eigen99_extraction"

#gwas1=fread("~/projects/MVA_output/02_PD_gwas/gwas/nallsEtAl2019_excluding23andMe.txt",data.table=F)
#gwas2=fread("~/projects/MVA_output/02_PD_gwas/gwas/IPDGC_AAO_GWAS_sumstats_april_2018.txt",data.table=F)

#res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)
#res.HDL
