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

	if (ncol(nsigs)==5){
		cs_qc$GIP1_avarage_Z2[i]=nsigs[ind_gip1,5]
		cs_qc$max_avarage_Z2[i]=max(nsigs[ind,5])
	}
	

	cs_qc$max_nsnps[i]=max(nsigs[ind,2])
	cs_qc$max_nsigs[i]=max(nsigs[ind,3])
	cs_qc$max_ratio[i]=max(nsigs[ind,4])
}


#colnames(cs_qc)
# [1] "N_cluster"        "n_traits"         "study_ids"        "traits"          
# [5] "EFOs"             "h2s"              "intercepts"       "alfa"            
# [9] "expected_h2_gip1" "expected_rgs"     "ChatGPT_name"     "intercept_min"   
#[13] "intercept_max"    "GIP1_intecept"    "GIP1_sig"         "GIP1_nsnps"      
#[17] "GIP1_ratio"       "GIP1_avarage_Z2"  "max_avarage_Z2"   "max_nsnps"       
#[21] "max_nsigs"        "max_ratio"  

#heritability

lst=strsplit(cs_qc$h2s,";")
max_h2=unlist(lapply(lst,max))

summary(cs_qc$GIP1_h2/as.numeric(max_h2))

#inflation
summary(cs_qc$GIP1_intecept/cs_qc$intercept_max)
t.test(cs_qc$GIP1_intecept/cs_qc$intercept_max-1)

# sig loci

summary(cs_qc$GIP1_nsnps/cs_qc$max_nsnps)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01201 0.59667 0.74183 0.70047 0.96996 0.99297

summary(lm(cs_qc$GIP1_sig~cs_qc$max_nsigs))
summary(lm(cs_qc$GIP1_ratio~cs_qc$max_ratio))

l=cs_qc$GIP1_sig*(cs_qc$max_nsnps/cs_qc$GIP1_nsnps)
summary(lm(l~cs_qc$max_nsigs))


l=cs_qc$GIP1_sig*(cs_qc$max_nsnps/cs_qc$GIP1_nsnps)*cs_qc$GIP1_intecept/cs_qc$intercept_max
summary(lm(l~cs_qc$max_nsigs))

summary(lm(cs_qc$GIP1_avarage_Z2~cs_qc$max_avarage_Z2))
summary(cs_qc$GIP1_avarage_Z2/cs_qc$max_avarage_Z2)