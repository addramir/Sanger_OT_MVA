setwd("Projects/Sanger_OT_MVA/06_PD_fm/")

library(data.table)
full=fread("nall_l2gs_v1.csv")

th=fread("nalls_l2g_predictions_270423.csv")

table(full$gene_id%in%th$gene_id)

ind=which(full$gene_id%in%th$gene_id)

full=full[ind,]
full=full[!(is.na(full$y_proba_full_model)),]

ind=match(full$gene_id,th$gene_id)
th=th[ind,]
table(th$gene_id==full$gene_id)

table(is.na(full$y_proba_full_model))

plot(y=full$y_proba_full_model,x=th$y_proba,xlab="top hits",ylab="pd gwas")
abline(b=1,a=0)

cor(full$y_proba_full_model,th$y_proba)
summary(lm(full$y_proba_full_model~th$y_proba))

table(full$y_proba_full_model>0.5)
table(th$y_proba>0.5)
