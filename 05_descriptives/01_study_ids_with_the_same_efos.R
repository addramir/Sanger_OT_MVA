setwd("~/Projects/Sanger_OT_MVA/05_descriptives/")

library(data.table)

si=fread("20230518_study_index.csv",data.table=F)
si_s=si[si$has_sumstats==TRUE,]

si_sne=si_s[si_s$trait_efos!="[]",]

nrow(si_s)-nrow(si_sne)
#945

l=table(si_sne$trait_efos)
length(l)
#[1] 2830
table(l==1)

#FALSE  TRUE 
#1124  1706 
table(l==2)

#FALSE  TRUE 
#2297   533 
table(l>2)

#FALSE  TRUE 
#2239   591 
table(l>1&l<=50)
#FALSE  TRUE 
#1714  1116 
table(l>30)
#FALSE  TRUE 
#2819    11 

names(l)[l>30]
#[1] "['EFO_0000546']" "['EFO_0004530']" "['EFO_0004725']" "['EFO_0007010']" "['EFO_0007814']"
#[6] "['EFO_0007874']" "['EFO_0007937']" "['EFO_0008111']" "['EFO_0010118']" "['EFO_0010226']"
#[11] "['EFO_0010228']"

problematic_efos=names(l)[l>30]
table(si_sne$trait_efos%in%problematic_efos)
#FALSE  TRUE 
#5914  2037 

l1=l[l>1 & l<=30]
length(l1)
#1113

table(si_sne$trait_efos%in%names(l1))
#FALSE  TRUE 
#3743  4208

hist(l1,n=30,xlim = c(1,30),main="Histogram of number of traits with the same EFO",xlab="Number of traits with the same EFO")
summary(as.numeric(l1))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   2.000   3.000   3.781   4.000  30.000 






