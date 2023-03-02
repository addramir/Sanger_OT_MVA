setwd("~/Projects/Sanger_OT_MVA/01_draft/")
library("data.table")
gcor=fread("gcors_13x13.csv",data.table=F)
rownames(gcor)=gcor[,1]
gcor=gcor[,-1]
gcor=as.matrix(gcor)
gcor[gcor>1]=1
library(corrplot)
corrplot(gcor,method = "square",tl.cex = 0.5)

gcor[gcor<0.9]=0
corrplot(gcor,method = "square",tl.cex = 0.5,order = "hclust")

gcor=fread("gcor.csv",data.table=F)
rownames(gcor)=gcor[,1]
gcor=gcor[,-1]
gcor=as.matrix(gcor)
corrplot(gcor,method = "square",tl.cex = 0.8)

phe=fread("phen_corr.csv",data.table=F)
phe=as.matrix(phe)
corrplot(phe,method = "square",tl.cex = 0.8)