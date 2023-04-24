library("data.table")
library("HDL")

setwd("~/projects/MVA_output/02_PD_gwas/gwas/")
lst=list.files()
n=length(lst)

LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_SVD_eigen99_extraction"

out_rg=array(NA,c(n,n))
out_pval=array(NA,c(n,n))
h2=rep(NA,n)
h2_se=rep(NA,n)
for (i in 1:(n-1)){
	for (j in (i+1):n){
		gwas1=fread(lst[i],data.table=F)
		gwas2=fread(lst[j],data.table=F)
		res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)

		out_rg[i,j]=out_rg[j,i]=res.HDL$rg
		out_pval[i,j]=out_pval[j,i]=res.HDL$P
		h2[i]=res.HDL$estimates.df[1,"Estimate"]
		h2[j]=res.HDL$estimates.df[2,"Estimate"]
		h2_se[i]=res.HDL$estimates.df[1,"se"]
		h2_se[j]=res.HDL$estimates.df[2,"se"]
	}
}

save(list=c("out_rg","out_pval","h2","h2_se"),file="~/projects/MVA_output/02_PD_gwas/v1_hm3_matrix.RData")

######

LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_hapmap2_SVD_eigen99_extraction"

out_rg=array(NA,c(n,n))
out_pval=array(NA,c(n,n))
h2=rep(NA,n)
h2_se=rep(NA,n)
for (i in 1:(n-1)){
	for (j in (i+1):n){
		gwas1=fread(lst[i],data.table=F)
		gwas2=fread(lst[j],data.table=F)
		res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)

		out_rg[i,j]=out_rg[j,i]=res.HDL$rg
		out_pval[i,j]=out_pval[j,i]=res.HDL$P
		h2[i]=res.HDL$estimates.df[1,"Estimate"]
		h2[j]=res.HDL$estimates.df[2,"Estimate"]
		h2_se[i]=res.HDL$estimates.df[1,"se"]
		h2_se[j]=res.HDL$estimates.df[2,"se"]
	}
}

save(list=c("out_rg","out_pval","h2","h2_se"),file="~/projects/MVA_output/02_PD_gwas/v1_hm3_matrix.RData")





gwas1=fread("~/projects/MVA_output/02_PD_gwas/gwas/NEALE2_20107_11.txt",data.table=F)
gwas2=fread("~/projects/MVA_output/02_PD_gwas/gwas/NEALE2_20110_11.txt",data.table=F)

res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)
res.HDL
