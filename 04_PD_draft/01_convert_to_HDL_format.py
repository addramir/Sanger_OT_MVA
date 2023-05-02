import numpy as np
import pandas as pd
import re
import os


os.chdir("/home/yt4/projects/Sanger_OT_MVA/04_PD_draft/")
import core_functions as CF

##### Global variables
path_out="~/projects/MVA_output/02_PD_gwas/gwas/"
max_ram_to_use="60g"

list_of_ids=["FINNGEN_R6_PDSTRICT"                           
,"FINNGEN_R6_PDSTRICT_EXMORE"                    
,"FINNGEN_R6_PD_DEMENTIA"                        
,"FINNGEN_R6_PD_DEMENTIA_EXMORE"                 
,"FINNGEN_R6_G6_PARKINSON"                       
,"FINNGEN_R6_G6_PARKINSON_EXMORE"                
,"FINNGEN_R6_G6_PARKINSON_INCLAVO"               
,"FINNGEN_R6_G6_PARKSCND"                        
,"NEALE2_20002_1262"                             
,"NEALE2_20107_11"                               
,"NEALE2_20110_11"                               
,"NEALE2_20111_11"                               
,"GCST90000014"                                  
,"GCST90000015"                                  
,"SAIGE_332"                                     
,"FINNGEN_R6_OTHER_DRUGADVERS_SECONDA_PARKINSONI"
,"FINNGEN_R6_PD2ND"                              
,"FINNGEN_R6_PD2ND_EXMORE"                       
,"SAIGE_966"]

CF.save_gwas_hm3_snps_HDL(list_of_ids=list_of_ids, max_ram_to_use=max_ram_to_use, path_out=path_out)
#CF.save_gwas_all_snps_HDL(list_of_ids=list_of_ids, max_ram_to_use=max_ram_to_use, path_out=path_out)



path_out="~/projects/MVA_output/02_PD_gwas/full_gwas/"
list_of_ids=["NEALE2_20002_1262",
"NEALE2_20110_11","SAIGE_332"]
CF.save_gwas_all_snps_HDL(list_of_ids=list_of_ids, max_ram_to_use=max_ram_to_use, path_out=path_out)

list_of_ids=["NEALE2_20107_11",
"FINNGEN_R6_G6_PARKINSON_INCLAVO"]
CF.save_gwas_all_snps_HDL(list_of_ids=list_of_ids, max_ram_to_use=max_ram_to_use, path_out=path_out)


path_out="~/projects/MVA_output/02_PD_gwas/full_gwas/"
list_of_ids=["NEALE2_20107_11","NEALE2_20110_11",
"FINNGEN_R6_G6_PARKINSON_INCLAVO","FINNGEN_R6_G6_PARKINSON"]
CF.save_gwas_all_snps_HDL(list_of_ids=list_of_ids, max_ram_to_use=max_ram_to_use, path_out=path_out)