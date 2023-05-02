library("data.table")

setwd("~/projects/MVA_output/02_PD_gwas/full_gwas/")
load("~/projects/MVA_output/02_PD_gwas/5_Traits_Z_matrix.RData")
Z=out

source("~/projects/Sanger_OT_MVA/04_PD_draft/shared_heredity.R")
library(Rcpp)
sourceCpp("~/projects/Sanger_OT_MVA/04_PD_draft/linear_combination_v4.cpp")


#[1] "nallsEtAl2019_excluding23andMe.txt"  "NEALE2_20107_11.txt"                
#[3] "NEALE2_20110_11.txt"                 "FINNGEN_R6_G6_PARKINSON_INCLAVO.txt"
#[5] "FINNGEN_R6_G6_PARKINSON.txt" 

ncases=c(33674,8043,5443,2562,2496)
nctotal=c(482730,312104,330077,260405,260405)
prev=ncases/nctotal
neff=4*prev*(1-prev)*nctotal

gcor=as.matrix(array(1,c(5,5)))

h2_set=0.0162
gcov=h2_set*(sqrt(neff%o%neff)/(max(neff)))

colnames(gcov)=rownames(gcov)=colnames(phe)=rownames(phe)=lst
sh=shared_heredity(CovGenTr = gcov, CorPhenTr = phe)

a=sh$GIPs$GIP_coeff[,1]

Z=as.data.frame(Z)
Z=na.omit(Z)
Zs=Z[,c(2,5,8,11,14)]
N=Z[,c(3,6,9,12,15)]
eaf=Z[,4]

result1 <- GWAS_linear_combination_Z_basedC(a=a, Z=as.matrix(Zs), covm=as.matrix(phe), N=as.matrix(N), eaf=eaf)

res=cbind(Z[,c(1,4)],result1[,1:3])

ref=fread("~/projects/MVA_output/02_PD_gwas/full_gwas/NEALE2_20107_11.txt",data.table=F)

table(res$SNP%in%ref$SNP)
ind=match(res$SNP,ref$SNP)
res[,"A1"]=ref[ind,"A1"]
res[,"A2"]=ref[ind,"A2"]
res[,"chr"]=ref[ind,"chr"]
res[,"pos"]=ref[ind,"pos"]
colnames(res)[2]="eaf"


z=res$b/res$se
pval=pchisq(z^2,df=1,low=F)
res=cbind(res,pval)

fwrite(x=res,file="~/projects/MVA_output/02_PD_gwas/full_gwas/GIP1_5_traits.txt",sep="\t")
