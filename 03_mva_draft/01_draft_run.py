import scipy as sc
from scipy import stats
import numpy as np
import pandas as pd
import re
import os


os.chdir("/home/yt4/projects/Sanger_OT_MVA/03_mva_draft/")
import core_functions as CF

##### Global variables
max_in_mva=10
path_to_save="/home/yt4/projects/MVA_output/01_clusters_draft/output/"
max_ram_to_use="40g"

##### Load all data
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

##### Main loop

clst=1
#for clst in range(1,number_of_clusters+1):
for clst in range(1,number_of_clusters+1):
    print(str(clst))
    subclst=clusters[clusters.iloc[:,1]==clst]    
    if len(subclst)>1:
        list_of_ids=list(subclst.iloc[:,0])
        phe=CF.phe_corr(list_of_ids,max_ram_to_use)
        list_for_mva=CF.list_of_mva_traits_based_on_phe(phe=phe,max_in_mva=max_in_mva,list_of_ids=list_of_ids,h2=h2)
        ind = [list_of_ids.index(elem) for elem in list_for_mva] 
        phe_mva=phe[np.ix_(ind,ind)]
        ind = [list_of_traits_in_gcor.index(elem) for elem in list_for_mva] 
        gcor_mva=np.array(GCOR.iloc[ind,ind])
        ind = [list_of_traits_in_h2.index(elem) for elem in list_for_mva] 
        h2_mva=np.array(h2.iloc[ind,1])
        Z,N,eaf,DF=CF.prepare_Z_N_eaf(list_for_mva,max_ram_to_use)
        MVA=CF.GIP1_lin_comb_Z_based(Z=Z,covm=phe_mva,eaf=eaf,N=N, gcor=gcor_mva, h2=h2_mva)
        df=pd.concat([DF,MVA], axis=1)       
        #saving data
        folder_path=path_to_save+str(clst)+"/"
        if not(os.path.exists(folder_path)):
            os.mkdir(folder_path)        
        df.to_csv(folder_path+"/GIP1.csv",index=False)
        np.savetxt(folder_path+"/gcor.csv", gcor_mva, delimiter=',')
        np.savetxt(folder_path+"/phe.csv", phe_mva, delimiter=',')
        np.savetxt(folder_path+"/h2.csv", h2_mva, delimiter=',')
        np.savetxt(folder_path+"/h2.csv", h2_mva, delimiter=',')
        dd=pd.DataFrame(list_for_mva,columns=["StudyID"])
        dd.to_csv(folder_path+"/list_of_traits.csv",index=False)

