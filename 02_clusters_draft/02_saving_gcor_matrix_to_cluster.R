setwd("~/projects/Sanger_OT_MVA/02_clusters_draft")
load("All_rg.Rdata")

df=Correlation_matrix

df[df>1]=1
df[df<(-1)]=-1
diag(df)=1
table(is.na(df))

l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
l2=apply(df,MAR=2,FUN=function(x){sum(is.na(x))})/nrow(df)
table(l1==l2)
df=df[l1<0.005,l1<0.005]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
hist(l1)
df=df[l1<0.001,l1<0.001]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
hist(l1)
df=df[l1<0.0001,l1<0.0001]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
table(is.na(df))
dim(df)

library(data.table)

fwrite(x=df,file="~/projects/MVA_output/01_clusters_draft/df_1222.csv",sep=",",quote=F,col.names=T,row.names=T)