setwd("~/projects/Sanger_OT_MVA/03_mva_draft/")
library(data.table)
cs_qc=fread("20230523_ckuster_summary_for_QC.txt",data.table=F)

path_to_save="/home/yt4/projects/MVA_output/01_clusters_draft/output/"
i=1
for (i in 1:nrow(cs_qc)){
	nsigs=fread(paste0(path_to_save,cs_qc[i,1],"/number_of_sigs.csv"),data.table=F)

	ind=which((nsigs[,1]!="GIP1"))
	ind_gip1=which((nsigs[,1]=="GIP1"))

	cs_qc$GIP1_sig[i]=nsigs[ind_gip1,3]
	cs_qc$GIP1_nsnps[i]=nsigs[ind_gip1,2]
	cs_qc$GIP1_ratio[i]=nsigs[ind_gip1,4]

	cs_qc$max_nsnps[i]=max(nsigs[ind,2])
	cs_qc$max_nsigs[i]=max(nsigs[ind,3])
	cs_qc$max_ratio[i]=max(nsigs[ind,4])
}



summary(cs_qc$GIP1_nsnps/cs_qc$max_nsnps)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01201 0.59667 0.74183 0.70047 0.96996 0.99297

summary(lm(cs_qc$GIP1_sig~cs_qc$max_nsigs))
summary(lm(cs_qc$GIP1_ratio~cs_qc$max_ratio))

l=cs_qc$GIP1_sig*(cs_qc$max_nsnps/cs_qc$GIP1_nsnps)
summary(lm(l~cs_qc$max_nsigs))


l=cs_qc$GIP1_sig*(cs_qc$max_nsnps/cs_qc$GIP1_nsnps)*cs_qc$GIP1_intecept/cs_qc$intercept_max
summary(lm(l~cs_qc$max_nsigs))