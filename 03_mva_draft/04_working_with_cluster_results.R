setwd("C://Projects/Sanger_OT_MVA/03_mva_draft/")
library(data.table)
x=fread("mva_results_with_clst_name.txt",data.table=F,sep="\t")

table(!grepl(pattern = ";",x$EFOs) & x$EFOs!="")

lst=strsplit(x$h2s,";")

i=1
out=NULL
for (i in 1:208){
  h2s=max(as.numeric(lst[[i]]))
  h2_exp=x$expected_h2_gip1[i]
  out=c(out,h2_exp/h2s)
}

summary(out)
