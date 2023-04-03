import numpy as np
from scipy.stats import chi2
import pandas as pd

def phe_corr(list_of_ids):
    from pyspark.sql import SparkSession
    import pyspark.sql.functions as f
    import pyspark.sql.types as t
    from pyspark.sql import Window
    import scipy as sc
    from scipy import stats
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2

    #
    global spark
    spark = (
        SparkSession.builder
        .master('local[*]')
        .config('spark.driver.memory', '15g')
        .appName('spark')
        .getOrCreate()
    )

    gw=list_of_ids[0]
    print("N "+str(0)+": "+gw)
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
    Z=gwas.select("id","zscore")

    for gw in list_of_ids[1:]:
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
    #print(corr_matrix_array)
    #np.linalg.det(corr_matrix_array)

    #np.savetxt('phen_corr.csv', corr_matrix_array, delimiter=',')
    return corr_matrix_array


def GIP1_lin_comb_Z_based(Z, covm, eaf, N, gcor, h2):
    import numpy as np
    from scipy.stats import chi2
    import pandas as pd

    gcov=np.sqrt(np.outer(h2, h2))*gcor

    a=np.linalg.eig(gcov)[1][:,0]
    vary = np.sum(covm * np.outer(a, a))
    a=a/np.sqrt(vary)
    vary = np.sum(covm * np.outer(a, a))
    if a[0]<0:
        a=-a
    print("GIP1 coeffs are:")
    print(a)    
    MVA=GWAS_linear_combination_Z_based(a=a,Z=Z,covm=covm,eaf=eaf,N=N)

    return MVA


def GWAS_linear_combination_Z_based(a, Z, covm, eaf, N):
    """
    This function performs a linear combination of genome-wide association studies (GWAS) based on their Z-scores,
    phenotypic correlations, and effective allele frequencies (EAF).

    Parameters:
    -----------
    a : numpy array, shape (1, m)
        The array of linear coefficients for the GWAS, where `m` is the number of GWAS.

    Z : numpy array, shape (M, m)
        The array of Z-scores for the GWAS, where `M` is the number of SNPs and `m` is the number of GWAS.

    covm : numpy array, shape (m, m)
        The phenotypic correlation matrix between the GWAS.

    eaf : numpy array, shape (1, M)
        The array of effective allele frequencies for each SNP, where `M` is the number of SNPs.

    N : numpy array, shape (M, m)
        The array of total sample sizes for each SNP in each GWAS, where `M` is the number of SNPs and `m` is the
        number of GWAS.

    Returns:
    --------
    out : pandas DataFrame
        The DataFrame containing the results of the linear combination of the GWAS. The columns are:

        - `beta` : numpy array, shape (M,)
            The array of effect sizes for each SNP.

        - `se` : numpy array, shape (M,)
            The array of standard errors for each effect size.

        - `pval` : numpy array, shape (M,)
            The array of p-values for each effect size.

        - `n_total` : float
            The effective sample size of the linear combination of the GWAS.
    """
    import numpy as np
    from scipy.stats import chi2
    import pandas as pd
    # Calculate standard errors from Z-scores
    SE = np.sqrt(1 / (Z ** 2 + N))
    # Calculate betas from Z-scores
    BETA = Z * SE
    # Calculate variance of y
    vary = np.sum(covm * np.outer(a, a))
    # Calculate linear combination of betas
    b = np.matmul(BETA, a)
    # Calculate variance of b
    def varb_lin(i):
        sei=SE[i,:]
        out=np.sum(covm*np.outer(sei,sei)*np.outer(a,a))
        return out
    varb=[varb_lin(i) for i in range(0,Z.shape[0])]
    varb=np.array(varb)
    seb=np.sqrt(varb)
    # Calculate effective sample size
    ind = np.where(np.abs(b / seb) < 2)
    N_eff = np.median(1 / (seb[ind] ** 2))
    # Calculate standard error and beta for variance normalization
    varg = 2 * (1 - eaf) * eaf
    b_norm = b / np.sqrt(varg)
    seb_norm = seb / np.sqrt(varg)
    # Calculate p-values using chi-squared distribution
    pval = chi2.sf((b_norm / seb_norm) ** 2, 1)
    # Create output DataFrame
    out = pd.DataFrame({'beta': b_norm, 'se': seb_norm, 'pval': pval, 'n_total': N_eff})
    return out


def GWAS_linear_combination_Z_based_OLD(a, Z, covm, eaf, N):
    """
    This fucntion performs linear combination of GWAS based on their Z-scores,phenotypic corrs, and eaf
    m is the number of GWAS, M is the number of SNPs
    a is the np array of linera coefficients 1xm
    Z is np array of Z scores Mxm
    Z is np array of N total Mxm
    covm is phenotypic coreralashion matrix mxm
    eaf is the np array of EAF for each SNP 1xM
    """

    from scipy.stats import chi2
    import numpy as np
    import pandas as pd
    #
    m=Z.shape[1]
    #
    def se_from_z(i):
        z=Z[:,i]
        n=N[:,i]
        se_x=np.sqrt(1/(z**2+n))
        return se_x
    SE=[se_from_z(i) for i in range(0,m)]
    SE=np.array(SE)
    SE=np.transpose(SE)
    #   
    def beta_from_z(i):
        z=Z[:,i]
        se=SE[:,i]
        beta=z*se
        return beta
    BETA=[beta_from_z(i) for i in range(0,m)]
    BETA=np.array(BETA)
    BETA=np.transpose(BETA)
    #
    vary=np.sum(covm*np.outer(a,a))
    b=np.matmul(BETA,a)
    #
    def varb_lin(i):
        sei=SE[i,:]
        out=np.sum(covm*np.outer(sei,sei)*np.outer(a,a))
        return out
    varb=[varb_lin(i) for i in range(0,Z.shape[0])]
    varb=np.array(varb)
    seb=np.sqrt(varb)
    #
    b=b/np.sqrt(vary)
    seb=seb/np.sqrt(vary)
    #
    ind = np.where(np.abs(b / seb) < 2)
    N_eff = np.median(1 / (seb[ind] ** 2))
    #
    varg = 2 * (1 - eaf) * eaf
    b=b/np.sqrt(varg)
    seb=seb/np.sqrt(varg)
    #
    from scipy.stats.distributions import chi2
    pval=chi2.sf((b/seb)**2,1)
    out=pd.DataFrame({'beta': b, 'se': seb, 'pval':pval, 'n_total':N_eff})
    return out

