setwd("Desktop/")
library(data.table)
x=fread("MVA_rg_comparison.csv",data.table=F)

head(x)



x$rg[x$rg>1]=1
x$rg[x$rg<(-1)]=-1

delta=abs(x$rg-x$expected_rgs)/abs(x$expected_rgs)


sum(table(x$MVA)==2)



summary(delta)

hist(delta,n=20)

table(delta>0.05)

y=x[delta>0.05,]

unique(y$MVA)


ll=matrix(c(0.0266    , 0.02999918,
       0.02999918, 0.039),nrow=2)





