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
max_ram_to_use="25g"


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
    nsigs=CF.number_of_sig_snps(list_of_ids,max_ram_to_use)
    v=["GIP1",len(df),sum(df["pval"]<=5e-8),sum(df["pval"]<=5e-8)/len(df)]
    sigs=pd.concat([nsigs, pd.DataFrame([v])], axis=0)
    sigs.columns=["ID","N","Nsig","Ratio"]
    sigs.to_csv(path_to_save+clst+"/number_of_sigs.csv",index=False)

