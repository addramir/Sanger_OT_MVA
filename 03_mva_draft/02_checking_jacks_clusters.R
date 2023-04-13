setwd("Desktop/")
load("All_rg.Rdata")
dim(Correlation_matrix)
df=Correlation_matrix

diag(df)=1
df[df>1]=1
df[df<(-1)]=-1
library(data.table)
x=fread("updated_grps_final_noNA.csv",header=T,data.table=F)

max_clsuters=max(x$grps)

i=1
out=NULL
for (i in 1:max_clsuters){
  ind=which(x$grps==i)
  sids=x[ind,1]
  dff=df[sids,sids]
  out=c(out,eigen(dff)$values[1]/sum(eigen(dff)$values))
}

summary(out)

table(x[,1]%in%colnames(df))

dff=df[x[,1],x[,1]]
table(is.na(dff))

l1=apply(dff,MAR=1,FUN=function(x) sum(is.na(x)))
summary(l1)
dff=dff[l1<10,l1<10]
l1=apply(dff,MAR=1,FUN=function(x) sum(is.na(x)))
summary(l1)
dff=dff[l1<5,l1<5]
l1=apply(dff,MAR=1,FUN=function(x) sum(is.na(x)))
summary(l1)
dff=dff[l1<1,l1<1]
l1=apply(dff,MAR=1,FUN=function(x) sum(is.na(x)))
summary(l1)
dim(dff)

df=dff
i=1
out=NULL
y=x[x[,1]%in%colnames(df),]
for (i in 1:max_clsuters){
  ind=which(y$grps==i)
  sids=y[ind,1]
  dff=df[sids,sids]
  if (length(sids)>0){
    out=c(out,eigen(dff)$values[1]/sum(eigen(dff)$values))
  } else{
    out=c(out,1)
  }
    
  
}
summary(out)


dim(df)
save(list=c("df"),file="df_1222.RData"))

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
out_m=array(NA,c(700,2))
for (i in 190:700){
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

grps=cutree(l,k = 617)
l_my=table(grps)
l_j=table(x$grps)
t.test(as.vector(l_j),as.vector(l_my))

summary(as.vector(l_j))
summary(as.vector(l_my))
table(as.vector(l_my)>5)
table(as.vector(l_j)>5)













