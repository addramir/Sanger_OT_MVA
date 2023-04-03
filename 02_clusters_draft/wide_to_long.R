setwd("../Desktop/ss/")
library(data.table)
x=fread("All_rg.txt",data.table=F)
x=x[,1:5]
x[,1]=gsub(x=x[,1],replacement = "",pattern = "munged_sumstats/")
x[,2]=gsub(x=x[,2],replacement = "",pattern = "munged_sumstats/")
x[,1]=gsub(x=x[,1],replacement = "",pattern = ".sumstats.gz")
x[,2]=gsub(x=x[,2],replacement = "",pattern = ".sumstats.gz")
head(x[,1:3])

p1=x[,1]
p2=x[,2]
x[,1]=p2
x[,2]=p1
# Pivot data into wide format
library(tidyr)
wide_df <- pivot_wider(x[,1:3], names_from = "p2", values_from = "rg", values_fill = NA)

wide_df=as.matrix(wide_df)

rnms=wide_df[,1]
wide_df=wide_df[,-1]
dim(wide_df)
rownames(wide_df)=rnms
length(unique(colnames(wide_df)))
length(unique(rownames(wide_df)))
rownames(wide_df)[!(rownames(wide_df)%in%colnames(wide_df))]
#"NEALE2_20116_1" "NEALE2_137"
wide_df=cbind(wide_df,V1=NA,V2=NA)
colnames(wide_df)[1248:1249]=rownames(wide_df)[!(rownames(wide_df)%in%colnames(wide_df))]
dim(wide_df)
wide_df=wide_df[rownames(wide_df),rownames(wide_df)]
wide_df[1:5,1:5]
df=as.matrix(apply(wide_df,MAR=1,FUN=as.numeric))
rownames(df)=colnames(df)
#
diag(df)=1
table(is.na(df))
i=1
out=NULL
for (i in (1:(ncol(df)-1))){
  for (j in ((i+1):ncol(df))){
    a1=df[i,j]
    a2=df[j,i]
    if(!is.na(a1)&!is.na(a2)&(a1!=a2)) print("mistake!")
    if (is.na(a1)) df[i,j]=a2
    if (is.na(a2)) df[j,i]=a1
    if(is.na(a1)&is.na(a2)) out=c(out,paste(colnames(df)[i],colnames(df)[j]))
  }
}
table(is.na(df))
write.table(out,file="list_of_NA.txt",sep=",",quote=F,col.names=F,row.names=F)
save(list=c("df"),file="matrix_v1.RData")



df1=df
df2=df
df1[upper.tri(df1)]=NA
df2[lower.tri(df2)]=NA
diag(df1)=NA
diag(df2)=NA

df2=t(df2)
df1[is.na(df1)]=df2[is.na(df1)]
df1[upper.tri(df1)]=t(df1[lower.tri(df1)])
diag(df1)=1
table(is.na(df1))
# 
l1=apply(df1,MAR=1,FUN=function(x){sum(is.na(x))})
l2=apply(df1,MAR=2,FUN=function(x){sum(is.na(x))})
