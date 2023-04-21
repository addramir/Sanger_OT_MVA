library("data.table")
library("HDL")

gwas1=fread("~/projects/MVA_output/02_PD_gwas/gwas/FINNGEN_R6_PD_DEMENTIA.txt",data.table=F)
gwas2=fread("~/projects/MVA_output/02_PD_gwas/gwas/FINNGEN_R6_PDSTRICT_EXMORE.txt",data.table=F)
LD.path <- "~/projects/MVA_output/02_PD_gwas/UKB_imputed_SVD_eigen99_extraction"
res.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 6)
res.HDL