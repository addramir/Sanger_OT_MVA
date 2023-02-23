from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np
import pandas as pd

#Spark initialization and configuration

global spark
spark = (
    SparkSession.builder
    .master('local[*]')
    .config('spark.driver.memory', '15g')
    .appName('spark')
    .getOrCreate()
)

#variant_annotation = spark.read.parquet("gs://genetics-portal-dev-data/22.09.0/outputs/lut/variant-index")
#variant_annotation = (variant_annotation
#               .withColumn("id", f.concat_ws("_", f.col("chr_id"), f.col("position"), f.col("ref_allele"), f.col("alt_allele")))
#               )

h2=pd.read_csv("h2.csv")
phe=pd.read_csv("phen_corr.csv",header=None)

SIDS=h2["study_id"].values

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

DF=DF.toPandas()

Z=DF[["z1","z2","z3","z4"]]
Z=np.array(Z)

N=DF[["n1","n2","n3","n4"]]
N=np.array(N)

a=[0.4980116,0.3785392,0.2401656,0.2578737]
a=np.array(a)

phe=np.array(phe)

eaf=DF["eaf"]
eaf=np.array(eaf)

from core_functions import GWAS_linear_combination_Z_based
#from core_functions import GWAS_linear_combination_Z_based_OLD

MVA=GWAS_linear_combination_Z_based(a=a,Z=Z,covm=phe,eaf=eaf,N=N)
#MVA2=GWAS_linear_combination_Z_based_OLD(a=a,Z=Z,covm=phe,eaf=eaf,N=N)

df=DF[["id","chrom","pos","ref","alt"]]
df=pd.concat([df,MVA], axis=1)

df.to_csv("/home/yt4/projects/SS_QC/CVD_MVA.csv",index=False)

