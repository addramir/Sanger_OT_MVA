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
#N1=t(apply(N,MAR=1,FUN=function(x){neff}))
eaf=Z[,4]

result1 <- GWAS_linear_combination_Z_basedC(a=a, Z=as.matrix(Zs), covm=as.matrix(phe), N=as.matrix(N), eaf=eaf)

res=cbind(Z[,c(1,4)],result1[,1:3])

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


#####
sss=c(
"rs356220",
"rs12752133",
"rs34311866",
"rs35603727",
"rs2696587")

res[match(sss,res$SNP),]
#               SNP    eaf           b          se        N   A1   A2  chr
#NA            <NA>     NA          NA          NA       NA <NA> <NA> <NA>
#2657348 rs12752133 0.0186  0.09722130 0.008509023 373922.3    T    C    1
#8009620 rs34311866 0.1958  0.03081664 0.002894268 373922.3    C    T    4
#8229990 rs35603727 0.0203  0.07584154 0.008143975 373922.3    A    G    1
#7581741  rs2696587 0.2173 -0.02750643 0.002784855 373922.3    G    T   17
#              pos         pval
#NA             NA           NA
#2657348 155205378 3.112406e-30
#8009620    951947 1.791630e-26
#8229990 156007988 1.247469e-20
#7581741  44222460 5.229875e-23
#rs356220     rs356220       t       c 0.3901 0.0000  0.3901  0.3901  0.2326
#rs12752133 rs12752133       t       c 0.0293 0.0114  0.0186  0.0415  0.6367
#rs34311866 rs34311866       t       c 0.8009 0.0051  0.7930  0.8042 -0.2085
#rs35603727 rs35603727       a       g 0.0272 0.0087  0.0203  0.0381  0.5264
#rs2696587   rs2696587       t       g 0.8074 0.0520  0.7827  0.9172  0.2284
#           StdErr   P-value Direction HetISq HetChiSq HetDf  HetPVal
#rs356220   0.0173 3.290e-41        +?    0.0    0.000     0 1.000000
#rs12752133 0.0510 8.220e-36        ++   89.2    9.294     1 0.002299
#rs34311866 0.0194 6.811e-27        --   55.1    2.225     1 0.135700
#rs35603727 0.0491 8.696e-27        ++    0.0    0.119     1 0.729800
#rs2696587  0.0230 3.741e-23        ++   42.7    1.745     1 0.186600
#           TotalSampleSize       pos chr
#rs356220            482730  90641340   4
#rs12752133          740998 155205378   1
#rs34311866          743135    951947   4
#rs35603727          743135 156007988   1
#rs2696587           743135  44222460  17

#               SNP    eaf           b          se        N   A1   A2  chr
#NA            <NA>     NA          NA          NA       NA <NA> <NA> <NA>
#2657348 rs12752133 0.0186  0.09618431 0.008298351 392204.9    T    C    1
#8009620 rs34311866 0.1958  0.03051788 0.002822014 392204.9    C    T    4
#8229990 rs35603727 0.0203  0.07402399 0.007940701 392204.9    A    G    1
#7581741  rs2696587 0.2173 -0.02750930 0.002715331 392204.9    G    T   17
#              pos         pval
#NA             NA           NA
#2657348 155205378 4.589667e-31
#8009620    951947 2.947957e-27
#8229990 156007988 1.140630e-20
#7581741  44222460 4.020808e-24

ind=c(1,4)
sh=shared_heredity(CovGenTr = gcov[ind,ind], CorPhenTr = phe[ind,ind])

a=sh$GIPs$GIP_coeff[,1]

result1 <- GWAS_linear_combination_Z_basedC(a=a, Z=as.matrix(Zs[,ind]), covm=as.matrix(phe[ind,ind]), N=as.matrix(N[,ind]), eaf=eaf)

res=cbind(Z[,c(1,4)],result1[,1:3])

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
res[match(sss,res$SNP),]


i=1
j=1
max_overlap=array(NA,c(length(neff),length(neff)))
for (i in 1:length(neff)){
	for (j in 1:length(neff)){
		n1=max(neff[i],neff[j])
		n2=min(neff[i],neff[j])
		max_overlap[i,j]=n2/sqrt(neff[j]*neff[i])
	}
}

ovpr=phe/max_overlap
ovpr[ovpr<0.9]=0




