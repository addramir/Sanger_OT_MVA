import scipy as sc
from scipy import stats
import numpy as np
import pandas as pd
import re
import os
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import squareform

os.chdir("/home/yt4/projects/Sanger_OT_MVA/03_mva_draft/")
import core_functions as CF


h2=pd.read_csv("~/projects/Sanger_OT_MVA/03_mva_draft/all_h2.csv",sep="\t")
h2["Trait"]=[re.sub(string=x,pattern=".results",repl="") for x in h2["Trait"]]
list_of_traits_in_h2=list(h2["Trait"])

GCOR=pd.read_csv("~/projects/MVA_output/01_clusters_draft/df_1222.csv",index_col=0)
list_of_traits_in_gcor=list(GCOR.columns)

#self checking
sum([elem in list_of_traits_in_h2 for elem in list_of_traits_in_gcor])

QC=pd.read_csv("/home/yt4/projects/Sanger_OT_MVA/03_mva_draft/QC_results.csv",sep=";")
QC=QC[QC["total_SNP"]>=2e6]
list_of_traits_in_qc=list(QC["study_id"])

clusters=pd.read_csv("~/projects/Sanger_OT_MVA/03_mva_draft/Study_and_cluster.csv")
list_of_traits_in_clusters=list(clusters.iloc[:,0])
sum([elem in list_of_traits_in_gcor for elem in list_of_traits_in_clusters])

clusters=clusters[clusters.iloc[:,0].isin(list_of_traits_in_gcor)]
clusters=clusters[clusters.iloc[:,0].isin(list_of_traits_in_qc)]

number_of_clusters=clusters["x"].max()

max_in_mva=3

clst=1
#for clst in range(1,number_of_clusters+1):
    subclst=clusters[clusters.iloc[:,1]==clst]
    
    if len(subclst)>1:
    
        list_of_ids=list(subclst.iloc[:,0])
    
        phe=CF.phe_corr(list_of_ids)
    
        if len(list_of_ids)>max_in_mva:
            dist = squareform(1 - np.power(phe, 2))
            l = linkage(dist, method='ward')
            grps = cut_tree(l, n_clusters=max_in_mva).flatten()
            out = []
            for i in range(0, max_in_mva):
                ids = list(np.array(list_of_ids)[grps == i])
                subset = h2[h2["Trait"].isin(ids)]
                si = subset.iloc[subset["Total Observed scale h2"].argmax(),0]
                out.append(si)
        else:
            out=list_of_ids

        list_for_mva=out
            
        ind = [list_of_ids.index(elem) for elem in list_for_mva] 
        phe_mva=phe[np.ix_(ind,ind)]

        ind = [list_of_traits_in_gcor.index(elem) for elem in list_for_mva] 
        gcor_mva=np.array(GCOR.iloc[ind,ind])

        ind = [list_of_traits_in_h2.index(elem) for elem in list_for_mva] 
        h2_mva=np.array(h2.iloc[ind,1])
        
        Z,N,eaf,DF=CF.prepare_Z_N_eaf(list_of_mva)




gw=SIDS[0]
print("N "+str(0)+": "+gw)
gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
gwas=spark.read.parquet(gw)
gwas=(gwas
    .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
    .withColumn("zscore", f.col("beta")/f.col("se"))
    )
gwas=(gwas
    .filter(f.col("eaf")>=0.001)
    .filter(f.col("eaf")<=0.999)
    )
DF=gwas.select("id","zscore","eaf","n_total",'chrom','pos','ref','alt')
DF=(DF
    .withColumnRenamed("zscore", "z1")
    .withColumnRenamed("n_total", "n1")
    )

i=2
for gw in SIDS[1:]:
    print("N: "+gw)
    gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
    gwas=spark.read.parquet(gw)
    gwas=(gwas
        .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
        .withColumn("zscore", f.col("beta")/f.col("se"))
        )
    gwas=(gwas
        .filter(f.col("eaf")>=0.001)
        .filter(f.col("eaf")<=0.999)
        )
    gwas=gwas.select("id","zscore","n_total")
    gwas=(gwas
        .withColumnRenamed("zscore", ("z"+str(i)))
        .withColumnRenamed("n_total", ("n"+str(i)))
    )
    i=i+1
    DF = DF.join(gwas, on = "id", how = "inner")


l=DF.join(variant_annotation, on="id", how="inner").distinct()
DF=l
DF=DF.toPandas()

sum(DF["z1"]**2>=29.71679)
sum(DF["z2"]**2>=29.71679)
sum(DF["z3"]**2>=29.71679)
sum(DF["z4"]**2>=29.71679)

Z=DF[["z1","z2","z3","z4"]]
Z=np.array(Z)

N=DF[["n1","n2","n3","n4"]]
N=np.array(N)

eaf=DF["eaf"]
eaf=np.array(eaf)

h2=np.array(h2["h2"])
gcor=np.array(gcor)
#a=[0.4980116,0.3785392,0.2401656,0.2578737]
#a=np.array(a)

phe=np.array(phe)




#from core_functions import GWAS_linear_combination_Z_based_OLD

#MVA=GWAS_linear_combination_Z_based(a=a,Z=Z,covm=phe,eaf=eaf,N=N)
#MVA2=GWAS_linear_combination_Z_based_OLD(a=a,Z=Z,covm=phe,eaf=eaf,N=N)
MVA=GIP1_lin_comb_Z_based(Z=Z,covm=phe,eaf=eaf,N=N, gcor=gcor, h2=h2)


df=DF[["id","chrom","pos","ref","alt","rs_id"]]
df=pd.concat([df,MVA], axis=1)

df.to_csv("/home/yt4/projects/SS_QC/CVD_MVA_GIP1_with_rsid.csv",index=False)

