import numpy as np
import pandas as pd
import re
import os


os.chdir("/home/yt4/projects/Sanger_OT_MVA/04_PD_draft/")
import core_functions as CF

##### Global variables
path_out="~/projects/MVA_output/02_PD_gwas/gwas/"
max_ram_to_use="60g"

list_of_ids=[ "FINNGEN_R6_PDSTRICT_EXMORE",
  "FINNGEN_R6_PD_DEMENTIA",
  "FINNGEN_R6_PD_DEMENTIA_EXMORE",
  "FINNGEN_R6_G6_PARKINSON",
  "FINNGEN_R6_G6_PARKINSON_EXMORE",
  "FINNGEN_R6_G6_PARKINSON_INCLAVO",
  "FINNGEN_R6_G6_PARKSCND",
  "NEALE2_20002_1262",
  "NEALE2_20107_11",
  "NEALE2_20110_11",
  "NEALE2_20111_11",
  "GCST90000014",
  "GCST90000015",
  "SAIGE_332"]


  