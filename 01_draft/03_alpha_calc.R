setwd("~/Projects/Sanger_OT_MVA/01_draft/")
source("shared_heredity.R")

gcor=read.table("gcor.csv",sep=",",header=T,row.names=1)
gcor=as.matrix(gcor)
phe=read.table("phen_corr.csv",sep=",",header=F)
phe=as.matrix(phe)
h2_table=read.table("h2.csv",sep=",",header=T)
h2=h2_table[,2]

gcov=gcor*sqrt(h2%o%h2)
colnames(gcov)=colnames(phe)=colnames(gcor)
rownames(gcov)=rownames(phe)=rownames(gcor)

L=shared_heredity(CovGenTr = gcov,CorPhenTr = as.matrix(phe))
L$alphas["OPTIM.max(Shared/Phen)",]

L=shared_heredity(CorGenTr = gcor,CorPhenTr = as.matrix(phe),h2 = h2)
L$alphas["OPTIM.max(Shared/Phen)",]
#L$GIPs$GIP_coeff[,"GIP1"]
L$res["OPTIM.max(Shared/Phen)","h2(Gen/Phen)"]
#L$GIPs$H2["GIP1"]
