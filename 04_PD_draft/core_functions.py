def save_gwas_hm3_snps_HDL(list_of_ids,max_ram_to_use,path_out):

    import pandas as pd
    from pyspark.sql import DataFrame, SparkSession
    import pyspark.sql.functions as f
    from pyspark.sql.types import *
    from pyspark.sql.window import Window
    import numpy as np

    global spark
    spark = (
        SparkSession.builder
        .master('local[*]')
        .config('spark.driver.memory', max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )


    variant_index=spark.read.parquet("gs://genetics-portal-dev-data/22.09.1/outputs/lut/variant-index")

    HM3_SNPs=spark.read.csv('gs://genetics-portal-dev-analysis/xg1/Configs/listHM3.txt')

    variant_index = (variant_index
               .withColumn("snpid", f.concat_ws("_", f.col("chr_id"), f.col("position"), f.col("ref_allele"), f.col("alt_allele")))
               )

    VI=variant_index.select(f.col("chr_id_b37"),f.col("position_b37"),f.col("rs_id"),f.col("snpid"))


    HM3_SNPs=HM3_SNPs.withColumnRenamed("_c0", "rs_id")
    HM3_SNPs=HM3_SNPs.join(VI,["rs_id"]).distinct()

    gw=list_of_ids[0]

    cts=['rs_id','ref','alt','beta','se','n_total']

    for i,gw in enumerate(list_of_ids):
        print("N "+str(i)+": "+gw)
        gwl='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
        gwas=spark.read.parquet(gwl)
        gwas=gwas.withColumn("snpid", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
        l=gwas.join(HM3_SNPs, on="snpid", how="inner").distinct()
        l=l.distinct()
        l=l.select(cts)
        L=l.toPandas()
        L.columns=['SNP', 'A2', 'A1', 'b', 'se', 'N']
        L.to_csv(path_out+gw+".txt",sep="\t",index=False)

    spark.stop()
    print("Done!")
