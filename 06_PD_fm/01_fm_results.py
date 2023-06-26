max_ram_to_use="10g"

from pyspark.sql import DataFrame, SparkSession
import pandas as pd
import pyspark.sql.functions as f
from pyspark.sql.types import *
from pyspark.sql.window import Window
import numpy as np
import os


def get_spark_session(max_ram_to_use):
    return (
        SparkSession.builder
        .master('local[*]')
        .config("spark.driver.memory", max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )


spark = get_spark_session(max_ram_to_use)



#coloc=spark.read.parquet("gs://genetics-portal-dev-staging/coloc/220408/coloc_processed_w_betas.parquet")
#coloc=spark.read.parquet("gs://genetics-portal-dev-analysis/dc16/nalls_et_al_coloc/coloc_processed.parquet")



#credset  finemap_snp  finemap_snp.csv  finemap_snp_filtered  top_loci  top_loci.json.gz

credset = spark.read.json("gs://genetics-portal-dev-analysis/yt4/PD_finemapping/results/credset")
credset=credset.toPandas()
credset.to_csv("credset.txt",sep="\t",index=False)

finemap_snp = spark.read.csv("gs://genetics-portal-dev-analysis/yt4/PD_finemapping/results/finemap_snp",header='true')
finemap_snp=finemap_snp.toPandas()
finemap_snp.to_csv("finemap_snp.txt",sep="\t",index=False)

finemap_snp_filtered = spark.read.csv("gs://genetics-portal-dev-analysis/yt4/PD_finemapping/results/finemap_snp_filtered",header='true')
finemap_snp_filtered=finemap_snp_filtered.toPandas()
finemap_snp_filtered.to_csv("finemap_snp_filtered.txt",sep="\t",index=False)

top_loci = spark.read.json("gs://genetics-portal-dev-analysis/yt4/PD_finemapping/results/top_loci")
top_loci=top_loci.toPandas()
top_loci.to_csv("top_loci.txt",sep="\t",index=False)