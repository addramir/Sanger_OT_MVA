library("data.table")
setwd("~/projects/MVA_output/02_PD_gwas/full_gwas/")
source("~/projects/Sanger_OT_MVA/04_PD_draft/shared_heredity.R")

#full
fl=fread("../PD_gwas/nallsEtAl2019_excluding23andMe_allVariants.tab",data.table=F)

l=strsplit(fl$SNP,split=":")

l=unlist(l)
chr=l[seq(1,length(l),by=2)]
chr=gsub(chr,pattern="chr",repl="")
chr=as.numeric(chr)

pos=l[seq(2,length(l),by=2)]
pos=as.numeric(pos)

maf=pmin(fl$freq,1-fl$freq)
fl=cbind(fl,chr,pos,maf)

fl_maf=fl[maf>=0.001,]

n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="p",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 32 12
n_f_full=cbind(n_f,trait="full")

#17510617

####
#with rsid
fl=fread("nallsEtAl2019_excluding23andMe.txt",data.table=F)
ref=fread("NEALE2_20107_11.txt",data.table=F)
ind=match(fl$SNP,ref$SNP)
table(fl$SNP==ref$SNP[ind])

fl=cbind(fl,pos=ref$pos[ind],chr=ref$chr[ind])
maf=pmin(fl$eaf,1-fl$eaf)
fl=fl[maf>=0.001,]
dim(fl)
#[1] 10782058       10
n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 27 10
n_f_full_rsid=cbind(n_f,trait="full_rsid")

####
#with rsid
fl=fread("GIP1_5_traits.txt",data.table=F)
head(fl)
maf=pmin(fl$eaf,1-fl$eaf)
fl=fl[maf>=0.001,]
dim(fl)
#[1] 8650201      10
n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 23 10

####
#with rsid
fl=fread("GIP1_5_traits_neff.txt",data.table=F)
head(fl)
maf=pmin(fl$eaf,1-fl$eaf)
fl=fl[maf>=0.001,]
dim(fl)
#[1] 8650201      10
n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 23 10


####
#with rsid
fl=fread("GIP1_4_traits.txt",data.table=F)
head(fl)
maf=pmin(fl$eaf,1-fl$eaf)
fl=fl[maf>=0.001,]
dim(fl)
#[1] 7992854       8

ind=match(fl$SNP,ref$SNP)
table(fl$SNP==ref$SNP[ind])
fl=cbind(fl,pos=ref$pos[ind],chr=ref$chr[ind])

n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 28 10

####
fl=fread("../generic-metal/METAANALYSIS1.TBL",data.table=F)

ind=match(fl$MarkerName,ref$SNP)
table(fl$MarkerName==ref$SNP[ind])

fl=cbind(fl,pos=ref$pos[ind],chr=ref$chr[ind])

maf=pmin(fl$Freq1,1-fl$Freq1)
fl=fl[maf>=0.001,]

fl=fl[!is.na(fl$chr),]
dim(fl)
#11760184       18
n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="P-value",pos="pos",snp="MarkerName", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 30 18           MarkerName Allele1 Allele2  Freq1 FreqSE MinFreq MaxFreq  Effect
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




n_f_ma=cbind(n_f,trait="MA")



n1=n_f_full_rsid[,c("SNP","pval","chr","pos","trait")]
n2=n_f_ma[,c("MarkerName","P-value","chr","pos","trait")]
colnames(n2)=colnames(n1)
N=rbind(n1,n2)
n_f=function_for_shlop_29_03_2020(locus_table=N,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait="trait")
dim(n_f)



####
fl=fread("~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_14.txt",data.table=F)

n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 29  9

####
fl=fread("~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_all.txt",data.table=F)

n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 24  9

####
fl=fread("~/projects/MVA_output/02_PD_gwas/full_gwas/20230518_MA_over_1-4.txt",data.table=F)

n_f=function_for_shlop_29_03_2020(locus_table=fl,p_value="pval",pos="pos",snp="SNP", delta=5e5,chr="chr",thr=5e-8,trait=NULL)
dim(n_f)
#[1] 24  9

