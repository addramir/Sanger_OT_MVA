setwd("Projects/Sanger_OT_MVA/06_PD_fm/")

library(data.table)
#full=fread("nall_l2gs_v1.csv")
full=fread("nalls_l2g_table.csv")
dim(full)
full=full[!(is.na(full$y_proba_full_model)),]

full_05=full[full$y_proba_full_model>=0.5]
length(unique(full_05$gene_id))
length(unique(full_05$pos))

length(unique(full$pos))

dup_ids=unique(full$gene_id[duplicated(full$gene_id)])

id=dup_ids[1]
for (id in dup_ids){
  ind=which(full$gene_id==id)
  ind_max=which.max(full$y_proba_full_model[ind])
  ind=ind[-ind_max]
  full=full[-ind,]
}

#####
th=fread("nalls_l2g_predictions_270423.csv",data.table=F)
length(unique(th$pos))
length(unique(th$gene_id))

th_05=th[th$y_proba>=0.5,]
length(unique(th_05$gene_id))
length(unique(th_05$pos))

dup_ids=unique(th$gene_id[duplicated(th$gene_id)])
id=dup_ids[1]
for (id in dup_ids){
  ind=which(th$gene_id==id)
  ind_max=which.max(th$y_proba[ind])
  ind=ind[-ind_max]
  th=th[-ind,]
}


table(full$gene_id%in%th$gene_id)

ind=which(full$gene_id%in%th$gene_id)

full=full[ind,]
length(unique(full$pos))

ind=match(full$gene_id,th$gene_id)
th=th[ind,]
table(th$gene_id==full$gene_id)

table(is.na(full$y_proba_full_model))

plot(y=full$y_proba_full_model,x=th$y_proba,xlab="top hits",ylab="pd gwas")
abline(b=1,a=0)

cor(full$y_proba_full_model,th$y_proba)
summary(lm(full$y_proba_full_model~th$y_proba))

table(full$y_proba_full_model>=0.5)
table(th$y_proba>=0.5)

ind1=which(full$y_proba_full_model>=0.5)
ind2=which(th$y_proba>=0.5)

length(intersect(ind1,ind2))
ind=unique(c(ind1,ind2))

plot(y=full$y_proba_full_model[ind],x=th$y_proba[ind],xlab="top hits",ylab="pd gwas")
abline(b=1,a=0)


