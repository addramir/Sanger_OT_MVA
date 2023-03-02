library(data.table)
x=fread("BETA_SE_for_IWV.csv",data.table=F)


IVW=function(B,SE,N){
	w=1/SE^2
	b=(B[,1]*w[,1]+B[,2]*w[,2])/(w[,1]+w[,2])
	se=sqrt(1/(w[,1]+w[,2]))
	n=N[,1]+N[,2]
	z=b/se
	pval=pchisq((z)^2,df=1,low=F)
	out=cbind(b,se,z,n,pval)
	out
}

l=IVW(B=x[,c("b1","b3")],SE=x[,c("se1","se3")],N=x[,c("n1","n3")])

gip1=fread("CVD_MVA_GIP1.csv",data.table=F)

ind=match(x[,"id"],gip1[,"id"])
table(x[,"id"]==gip1[ind,"id"])

z_gip1=gip1[,"beta"]/gip1[,"se"]

z_gip1=z_gip1[ind]

out=cbind(l[,"z"],z_gip1)

png(file="Z-Z_plot.png",width=10,height=10,units="in",res=300)
plot(x=l[,"z"],y=z_gip1,xlab="MA",ylab="GIP1")
dev.off()