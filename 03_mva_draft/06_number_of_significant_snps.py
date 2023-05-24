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
max_ram_to_use="100g"


def get_directories(directory_path):
    directories = []
    for entry in os.scandir(directory_path):
        if entry.is_dir():
            directories.append(entry.name)
    return directories

clusters=get_directories(path_to_save)
clst=clusters[0]

for clst in clusters:
    print(clst)
    df=pd.read_csv(path_to_save+clst+"/GIP1.csv")
    list_of_ids=pd.read_csv(path_to_save+clst+"/list_of_traits.csv")
    list_of_ids=list(list_of_ids["StudyID"])
    if (sum(df["pval"]<=5e-8)==0):
        nsigs=CF.number_of_sig_snps(list_of_ids,max_ram_to_use)
        v=["GIP1",len(df),sum(df["pval"]<=5e-8),sum(df["pval"]<=5e-8)/len(df)]
        sigs=pd.concat([nsigs, pd.DataFrame([v])], axis=0)
        sigs.columns=["ID","N","Nsig","Ratio"]
        sigs.to_csv(path_to_save+clst+"/number_of_sigs.csv",index=False)
    else:
        dfs=df[df["pval"]<=5e-8]
        dfs['id']= dfs.apply(lambda row: '_'.join([str(row['chrom']), str(row['pos']), row['ref'], row['alt']]), axis=1)
        list_of_snps=list(dfs["id"])
        nsigs=CF.number_of_sig_snps_and_mean_chi2(list_of_ids,max_ram_to_use,list_of_snps)
        v=["GIP1",len(df),sum(df["pval"]<=5e-8),sum(df["pval"]<=5e-8)/len(df),((dfs["beta"]/dfs["se"])**2).mean()]
        v=pd.DataFrame([v])
        v.columns=["ID","N","Nsig","Ratio","Avarage_Z2"]
        sigs=pd.concat([nsigs, v], axis=0)
        sigs.columns=["ID","N","Nsig","Ratio","Avarage_Z2"]
        sigs.to_csv(path_to_save+clst+"/number_of_sigs.csv",index=False)
        
