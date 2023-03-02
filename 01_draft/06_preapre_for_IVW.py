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


h2=pd.read_csv("h2.csv")
phe=pd.read_csv("phen_corr.csv",header=None)
gcor=pd.read_csv("gcor.csv",index_col=0)

SIDS=h2["study_id"].values

gw=SIDS[0]
print("N "+str(0)+": "+gw)
gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
gwas=spark.read.parquet(gw)
gwas=(gwas
    .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
    )
gwas=(gwas
    .filter(f.col("eaf")>=0.001)
    .filter(f.col("eaf")<=0.999)
    )
DF=gwas.select("id","beta","se","eaf","n_total",'chrom','pos','ref','alt')
DF=(DF
    .withColumnRenamed("beta", "b1")
    .withColumnRenamed("se", "se1")
    .withColumnRenamed("n_total", "n1")
    )

i=2
for gw in SIDS[1:]:
    print("N: "+gw)
    gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
    gwas=spark.read.parquet(gw)
    gwas=(gwas
        .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
        )
    gwas=(gwas
        .filter(f.col("eaf")>=0.001)
        .filter(f.col("eaf")<=0.999)
        )
    gwas=gwas.select("id","beta","se","n_total")
    gwas=(gwas
        .withColumnRenamed("beta", ("b"+str(i)))
        .withColumnRenamed("se", ("se"+str(i)))
        .withColumnRenamed("n_total", ("n"+str(i)))
    )
    i=i+1
    DF = DF.join(gwas, on = "id", how = "inner")

DF=DF.toPandas()

DF.to_csv("/home/yt4/projects/SS_QC/BETA_SE_for_IWV.csv",index=False)