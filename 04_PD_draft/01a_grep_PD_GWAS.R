setwd("~/Projects/Sanger_OT_GoldStandards/01_orig_GS/")

library(data.table)


si=fread("C://Projects/Sanger_OT_MVA/03_mva_draft/study_index.csv",data.table=F)

index=paste(si$trait_reported,si$n_initial,si$num_assoc_loci,sep = "_")
table(duplicated(index))
ind=which(duplicated(index) & si[,"num_assoc_loci"]>0)

dups=unique(index[ind])

y=si[index%in%dups,]

si_s=si[si$has_sumstats==TRUE,]
ind=grep(si_s$trait_reported,pattern = "Parkinson")
ind1=grep(si_s$trait_reported,pattern = "parkinson")

ind=unique(c(ind,ind1))

y=si_s[ind,]
