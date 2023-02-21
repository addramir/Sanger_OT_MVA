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
    .filter((f.col("zscore")**2)<=4)
    .filter(f.col("eaf")>=0.001)
    .filter(f.col("eaf")<=0.999)
    )
Z=gwas.select("id","zscore")

for gw in SIDS[1:]:
    print("N: "+gw)
    gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
    gwas=spark.read.parquet(gw)
    gwas=(gwas
        .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
        .withColumn("zscore", f.col("beta")/f.col("se"))
        )
    gwas=(gwas
        .filter((f.col("zscore")**2)<=4)
        .filter(f.col("eaf")>=0.01)
        .filter(f.col("eaf")<=0.99)
        )
    gwas=gwas.select("id","zscore")
    Z = Z.join(gwas, on = "id", how = "inner")


from pyspark.sql.functions import monotonically_increasing_id

# add unique identifier to column names
df = Z.toDF(*[f"{col}_{i}" for i, col in enumerate(Z.columns)])
# rename columns
df = df \
    .withColumnRenamed("id_0", "id") \
    .withColumnRenamed("zscore_1", "z1") \
    .withColumnRenamed("zscore_2", "z2") \
    .withColumnRenamed("zscore_3", "z3") \
    .withColumnRenamed("zscore_4", "z4")

# drop unique identifier column
df = df.drop("id_1")
Z=df.select('z1', 'z2', 'z3', 'z4')

from pyspark.ml.stat import Correlation
from pyspark.ml.feature import VectorAssembler

# Create a VectorAssembler to combine the columns into a single vector column
assembler = VectorAssembler(inputCols=['z1', 'z2', 'z3', 'z4'], outputCol='features')

# Apply the VectorAssembler to the DataFrame
data = assembler.transform(Z).select('features')

# Calculate the correlation matrix using the Pearson method
corr_matrix = Correlation.corr(data, 'features', method='pearson').head()
corr_matrix_array = corr_matrix[0].toArray()

# Print the correlation matrix
print(corr_matrix_array)
np.linalg.det(corr_matrix_array)

np.savetxt('phen_corr.csv', corr_matrix_array, delimiter=',')





