library("data.table")

setwd("~/projects/MVA_output/02_PD_gwas/full_gwas/")
load("~/projects/MVA_output/02_PD_gwas/5_Traits_Z_matrix.RData")
Z=out

source("~/projects/Sanger_OT_MVA/04_PD_draft/shared_heredity.R")
library(Rcpp)
sourceCpp("~/projects/Sanger_OT_MVA/04_PD_draft/linear_combination_v4.cpp")

ref=fread("~/projects/MVA_output/02_PD_gwas/full_gwas/NEALE2_20107_11.txt",data.table=F)

#[1] "nallsEtAl2019_excluding23andMe.txt"  "NEALE2_20107_11.txt"                
#[3] "NEALE2_20110_11.txt"                 "FINNGEN_R6_G6_PARKINSON_INCLAVO.txt"
#[5] "FINNGEN_R6_G6_PARKINSON.txt" 

ncases=c(33674,8043,5443,2562,2496)
nctotal=c(482730,312104,330077,260405,260405)
prev=ncases/nctotal


nqt_est=function(v,N){
  z=dnorm(x=qnorm(v,mean=0,sd=1,low=F),mean=0,sd=1)
  i=z/v
  Nqt=N*i^2*v/(1-v)
  Nqt
}

neff_est=function(v,N){
  neff=4*v*(1-v)*N
  neff
}

nqt_est_K=function(v,N,K){
  z=dnorm(x=qnorm(K,mean=0,sd=1,low=F),mean=0,sd=1)
  i=z/K
  Nqt=N*i^2*v*(1-v)/(1-K)^2
  Nqt
}


neff_1=neff_est(prev,nctotal)
neff_2=nqt_est(prev,nctotal)
neff_3=nqt_est_K(prev,nctotal,K=1/38)

Z=as.data.frame(Z)
Z=na.omit(Z)
Zs=Z[,c(2,5,8,11,14)]
N=Z[,c(3,6,9,12,15)]
#N1=t(apply(N,MAR=1,FUN=function(x){neff}))
eaf=Z[,c(4,7,10,13,16)]

VARG=2*eaf*(1-eaf)

SE=BETA=Zs

for (i in 1:ncol(Zs)){
  SE[,i]=sqrt(1/(VARG[,i]*N[,i]*prev[i]*(1-prev[i])))
  BETA[,i]=Zs[,i]*SE[,i]
}

function_for_MA=function(BETA,SE,phe){
  W=1/SE^2
  sW=sqrt(W)
  b_ma=se_ma=rep(NA,nrow(BETA))
  
  for (i in 1:nrow(BETA)){
      w=as.numeric(W[i,])
      se=as.numeric(SE[i,])
      a=w/sum(w)
      b_ma[i]=sum(BETA[i,]*a)
      se_ma[i]=sqrt(sum((se%o%se)*phe*(a%o%a)))
  }
  #b_ma <- rowSums(BETA * W) / rowSums(W)
  #se_ma=apply(sW,MAR=1,FUN=function(w) sqrt(1/sum((w%o%w)*phe)))


  pval=pchisq((b_ma/se_ma)^2,df=1,low=F)

  out=cbind(b_ma,se_ma,pval)
  out
}


L=function_for_MA(BETA[,c(1,4)],SE[,c(1,4)],phe[c(1,4),c(1,4)])
median(qchisq(L[,"pval"],df=1,low=F))/qchisq(0.5,df=1)

res=cbind(Z[,c(1,4)],L)

table(res$SNP%in%ref$SNP)
ind=match(res$SNP,ref$SNP)
res[,"A1"]=ref[ind,"A1"]
res[,"A2"]=ref[ind,"A2"]
res[,"chr"]=ref[ind,"chr"]
res[,"pos"]=ref[ind,"pos"]
colnames(res)[2]="eaf"

fwrite(x=res,file="~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_14.txt",sep="\t")



L=function_for_MA(BETA,SE,phe)
median(qchisq(L[,"pval"],df=1,low=F))/qchisq(0.5,df=1)
res=cbind(Z[,c(1,4)],L)

table(res$SNP%in%ref$SNP)
ind=match(res$SNP,ref$SNP)
res[,"A1"]=ref[ind,"A1"]
res[,"A2"]=ref[ind,"A2"]
res[,"chr"]=ref[ind,"chr"]
res[,"pos"]=ref[ind,"pos"]
colnames(res)[2]="eaf"
fwrite(x=res,file="~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_all.txt",sep="\t")



L=function_for_MA(BETA[,c(1,2,3,4)],SE[,c(1,2,3,4)],phe[c(1,2,3,4),c(1,2,3,4)])
median(qchisq(L[,"pval"],df=1,low=F))/qchisq(0.5,df=1)
res=cbind(Z[,c(1,4)],L)

table(res$SNP%in%ref$SNP)
ind=match(res$SNP,ref$SNP)
res[,"A1"]=ref[ind,"A1"]
res[,"A2"]=ref[ind,"A2"]
res[,"chr"]=ref[ind,"chr"]
res[,"pos"]=ref[ind,"pos"]
colnames(res)[2]="eaf"
fwrite(x=res,file="~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_1-4.txt",sep="\t")




sss=c(
"rs356220",
"rs12752133",
"rs34311866",
"rs35603727",
"rs2696587")
ind=match(sss,Z$SNP)


