library(qqman)
library(data.table)
x=fread("nallsEtAl2019_excluding23andMe.txt",data.table=F)
head(x)
dim(x)


ref=fread("NEALE2_20107_11.txt",data.table=F)
head(ref)
dim(ref)

table(x$SNP%in%ref$SNP)
ind=match(x$SNP,ref$SNP)
length(ind)
table(x$SNP==ref$SNP[ind])
x=cbind(x,chr=ref$chr[ind],pos=ref$pos[ind])

png(file="nalle_manhattan.png",height=500,width=1000,units="px")
x[,"chr"]=as.numeric(x[,"chr"])
manhattan(x, chr="chr", bp="pos", snp="SNP", p="pval")
dev.off()


x=fread("MA_Nalle_FinnGen.tsv",data.table=F)
head(x)
dim(x)
table(x$MarkerName%in%ref$SNP)
ind=match(x$MarkerName,ref$SNP)
table(x$MarkerName==ref$SNP[ind])

x=cbind(x,chr=ref$chr[ind],pos=ref$pos[ind])
x[,"chr"]=as.numeric(x[,"chr"])

y=na.omit(x)
dim(y)


png(file="MA_nalle_finngen_manhattan.png",height=500,width=1000,units="px")
manhattan(y, chr="chr", bp="pos", snp="MarkerName", p="P-value")
dev.off()