setwd("../Desktop/ss/")
load("matrix_v1.RData")

l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
l2=apply(df,MAR=2,FUN=function(x){sum(is.na(x))})/nrow(df)
table(l1==l2)
table(l1>0.9)

df=df[l1<0.5,l1<0.5]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
hist(l1)


df=df[l1<0.005,l1<0.005]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
hist(l1)

df=df[l1<0.001,l1<0.001]
l1=apply(df,MAR=1,FUN=function(x){sum(is.na(x))})/nrow(df)
hist(l1)


load("../Dropbox/ss/All_rg.Rdata")
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
hist(l1)
dim(df)


dst=as.dist(2*(1-abs(df)))
dst=as.dist((1-(df)^2))

l=hclust(d=dst,method="ward.D2")
df=df[l$order,l$order]
library(corrplot)
corrplot(df,method="square",tl.pos='n')

ll=eigen(df)
sum(ll$values[1:100])/sum(ll$values)
grps=cutree(l,k = 88)
table(grps)

clstr=12
names(grps)[grps==clstr]
corrplot(df[names(grps)[grps==clstr],names(grps)[grps==clstr]],method="square")
            
plot(cmdscale(dst,k = 2))


dff=df
dff[abs(dff)<0.7]=0
#dst=as.dist(2*(1-abs(df)))
dst=as.dist((1-(dff)^2))

l=hclust(d=dst,method="ward.D2")
dff=dff[l$order,l$order]


ll=eigen(dff)
sum(ll$values[1:190])/sum(ll$values)
grps=cutree(l,k = 500)
lll=table(grps)

i=190
out_m=array(NA,c(1000,2))
for (i in 190:1000){
  grps=cutree(l,k = i)
  lll=table(grps)  
  out=NULL
  for (j in 1:i){
    clstr=j
    dfff=df[names(grps)[grps==clstr],names(grps)[grps==clstr]]
    out=c(out,eigen(dfff)$values[1]/sum(eigen(dfff)$values))
  }
  out_m[i-189,1]=i
  out_m[i-189,2]=min(out)
}

plot(out_m[,2])

ind=which(out_m[,2]>=0.95)
out_m[min(ind),]

lll=c("FINNGEN_R6_PDSTRICT",
  "FINNGEN_R6_PDSTRICT_EXMORE",
  "FINNGEN_R6_PD_DEMENTIA",
  "FINNGEN_R6_PD_DEMENTIA_EXMORE",
  "FINNGEN_R6_G6_PARKINSON",
  "FINNGEN_R6_G6_PARKINSON_EXMORE",
  "FINNGEN_R6_G6_PARKINSON_INCLAVO",
  "FINNGEN_R6_G6_PARKSCND",
  "NEALE2_20002_1262",
  "NEALE2_20107_11",
  "NEALE2_20110_11",
  "NEALE2_20111_11",
  "GCST90000014",
  "GCST90000015",
  "SAIGE_332")

table(lll%in%names(grps))




grps=cutree(l,k = 619)
lll=table(grps) 
clstr=300
names(grps)[grps==clstr]
dfff=df[names(grps)[grps==clstr],names(grps)[grps==clstr]]
corrplot(dfff,method="square")
eigen(dfff)$values[1]/sum(eigen(dfff)$values)

table(lll>1)

hist(lll[lll>1],n=20)
table(lll)

plot(cmdscale(dst,k = 2))







