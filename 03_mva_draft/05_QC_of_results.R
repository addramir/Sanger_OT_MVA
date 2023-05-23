setwd("~/Projects/Sanger_OT_MVA/03_mva_draft/")
library(data.table)
rgs=fread("QC_results_rg_h2_formated.txt",data.table=F)
cluster_summary=fread("mva_results_with_clst_name.txt",data.table=F,sep="\t")


rgs$rg[rgs$rg>1]=1
rgs$rg[rgs$rg<(-1)]=-1

delta_rg=abs(rgs$rg-rgs$expected_rgs)/abs(rgs$expected_rgs)
delta_h2=abs(rgs$expected_h2-rgs$Total.Observed.scale.h2)/abs(rgs$expected_h2)
intercept=rgs$Intercept

ind_to_excl=which(delta_rg>0.1 | delta_h2>0.4 | intercept>1.3 | intercept<0.9)

clusters=unique(rgs$MVA)
clusters_to_exclude=unique(rgs$MVA[ind_to_excl])

clusters_passed_qc=clusters[!clusters%in%clusters_to_exclude]

ind=match(clusters_passed_qc,rgs$MVA)
table(rgs$MVA[ind]==clusters_passed_qc)
intercepts=rgs$Intercept[ind]

intercept_min=unlist(lapply(strsplit(cluster_summary$intercepts,split=";"),FUN=min))
intercept_max=unlist(lapply(strsplit(cluster_summary$intercepts,split=";"),FUN=max))
cluster_summary=cbind(cluster_summary,intercept_min,intercept_max)

clusters_passed_qc_s=paste0("./",clusters_passed_qc)
ind=match(cluster_summary$N_cluster,clusters_passed_qc_s)

table(clusters_passed_qc_s[ind]==cluster_summary$N_cluster)

cluster_summary$GIP1_intecept=NA
cluster_summary$GIP1_intecept=intercepts[ind]

CS_nona=na.omit(cluster_summary)
summary(CS_nona$GIP1_intecept/as.numeric(CS_nona$intercept_max))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7560  0.9627  1.0108  0.9896  1.0315  1.0650

fwrite(x=CS_nona,file="20230523_ckuster_summary_for_QC.txt",sep="\t")
