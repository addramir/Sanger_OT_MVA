setwd("~/projects/MVA_output/01_clusters_draft/output")

clstrs=list.dirs()
clstrs=clstrs[-1]

cl=clstrs[1]
#[1] "gcor.csv"           "GIP1.csv"           "h2.csv"
#[4] "list_of_traits.csv" "phe.csv"

library(data.table)
source("~/projects/Sanger_OT_MVA/03_mva_draft/shared_heredity.R")
si=fread("~/projects/Sanger_OT_MVA/03_mva_draft/study_index.csv",sep="\t",data.table=F)
h2_int=fread("~/projects/Sanger_OT_MVA/03_mva_draft/all_h2_intercepts.tsv",sep="\t",data.table=F)


out=as.data.frame(array(NA,c(length(clstrs),10)))
colnames(out)=c("N_cluster","n_traits","study_ids","traits","EFOs","h2s","intercepts","alfa","expected_h2_gip1","expected_rgs")

for (cl in clstrs){

	ltr=read.table(paste0(cl,"/list_of_traits.csv"),header=T)[,1]
	gcor=read.table(paste0(cl,"/gcor.csv"),header=F,sep=",")
	phe=read.table(paste0(cl,"/phe.csv"),header=F,sep=",")
	h2=read.table(paste0(cl,"/h2.csv"),header=F,sep=",")[,1]

	phe=as.matrix(phe)
	gcor=as.matrix(gcor)
	gcov=sqrt(h2%o%h2)*gcor

	gips=add_gips(gcovm=gcov,phem=phe)

	sub_si=si[match(ltr,si$study_id),]
	sub_h2=h2_int[match(ltr,h2_int$Trait),]

	efos=unique(sub_si[,"trait_efos"])
	efos=gsub(efos,pattern="'",replacement="",fixed=T)
	efos=gsub(efos,pattern="[",replacement="",fixed=T)
	efos=gsub(efos,pattern="]",replacement="",fixed=T)

	i=match(cl,clstrs)
	out[i,"N_cluster"]=cl
	out[i,"n_traits"]=length(ltr)
	out[i,"study_ids"]=paste(ltr,collapse=";")
	out[i,"traits"]=paste(sub_si[,"trait_reported"],collapse=";")
	out[i,"EFOs"]=paste(efos,collapse=";")
	out[i,"h2s"]=paste(h2,collapse=";")
	out[i,"intercepts"]=paste(sub_h2[,"Intercept"],collapse=";")
	out[i,"alfa"]=paste(gips$GIP_coeff[,1],collapse=";")
	out[i,"expected_h2_gip1"]=gips$H2["GIP1"]
	out[i,"expected_rgs"]=paste(gips$cor_g[1:length(ltr),"GIP1"],collapse=";")
}
fwrite(x=out,file="~/projects/Sanger_OT_MVA/03_mva_draft/mva_results_wo_clst_name.tsv",sep="\t")