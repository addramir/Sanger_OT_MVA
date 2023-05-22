
def list_of_mva_traits_based_on_phe(phe,max_in_mva,list_of_ids,h2):
    from scipy.cluster.hierarchy import linkage, cut_tree
    from scipy.spatial.distance import squareform
    import numpy as np
    from scipy.stats import chi2
    import pandas as pd
    
    if len(list_of_ids)>max_in_mva:
        dist = squareform(1 - np.power(phe, 2))
        l = linkage(dist, method='ward')
        grps = cut_tree(l, n_clusters=max_in_mva).flatten()
        out = []
        for i in range(0, max_in_mva):
            ids = list(np.array(list_of_ids)[grps == i])
            subset = h2[h2["Trait"].isin(ids)]
            si = subset.iloc[subset["Total Observed scale h2"].argmax(),0]
            out.append(si)
    else:
        out=list_of_ids

    return out

def prepare_Z_N_eaf(list_of_ids,max_ram_to_use):
    from pyspark.sql import SparkSession
    import pyspark.sql.functions as f
    import pyspark.sql.types as t
    from pyspark.sql import Window
    import scipy as sc
    from scipy import stats
    import numpy as np
    import pandas as pd
    from scipy.stats import chi2


    #Spark initialization and configuration

    global spark
    spark = (
        SparkSession.builder
        .master('local[*]')
        .config('spark.driver.memory', max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )

    variant_annotation = spark.read.parquet("gs://genetics-portal-dev-data/22.09.1/outputs/lut/variant-index")
    variant_annotation = (variant_annotation
                   .withColumn("id", f.concat_ws("_", f.col("chr_id"), f.col("position"), f.col("ref_allele"), f.col("alt_allele")))
                   )
    variant_annotation=variant_annotation.select("rs_id","id")

    gw=list_of_ids[0]
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
    DF=DF.join(variant_annotation, on="id", how="inner").distinct()
    DF=DF.toPandas()
    variant_annotation=0

    i=2
    for gw in list_of_ids[1:]:
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
        gwas=gwas.toPandas()
        DF=pd.merge(DF,gwas,on="id",how="inner")
        #DF = DF.join(gwas, on = "id", how = "inner")

    #DF=l.toPandas()

    spark.stop()

    ntraits=len(list_of_ids)
    Z=DF[["z"+str(elem) for elem in range(1,ntraits+1)]]
    Z=np.array(Z)

    N=DF[["n"+str(elem) for elem in range(1,ntraits+1)]]
    N=np.array(N)

    eaf=DF["eaf"]
    eaf=np.array(eaf)

    DF=DF[['rs_id','chrom', 'pos', 'ref', 'alt','eaf']]

    return Z,N,eaf,DF

def phe_corr(list_of_ids,max_ram_to_use):
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
        .config('spark.driver.memory', max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )

    gw=list_of_ids[0]
    print("N "+str(0)+": "+gw)
    gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
    gwas=spark.read.parquet(gw)
    #print(gwas.count())
    gwas=(gwas
        .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
        .withColumn("zscore", f.col("beta")/f.col("se"))
        )
    gwas=(gwas
        #.filter((f.col("zscore")**2)<=4)
        .filter(f.col("eaf")>=0.01)
        .filter(f.col("eaf")<=0.99)
        )
    Z=gwas.select("id","zscore")
    Z=Z.toPandas()
    Z.columns=["id","zscore1"]


    j=2
    for gw in list_of_ids[1:]:
        print("N: "+gw)
        gw='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
        gwas=spark.read.parquet(gw)
        #l=gwas.count()
        #if l>2e6:
        gwas=(gwas
            .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
            .withColumn("zscore", f.col("beta")/f.col("se"))
            )
        gwas=(gwas
            #.filter((f.col("zscore")**2)<=4)
            .filter(f.col("eaf")>=0.01)
            .filter(f.col("eaf")<=0.99)
            )
        gwas=gwas.select("id","zscore")
        #    
        gwas=gwas.toPandas()
        gwas.columns=["id","zscore"+str(j)]
        j=j+1
        #Z = Z.join(gwas, on = "id", how = "inner")
        Z=pd.merge(Z,gwas,on="id",how="inner")

    #from pyspark.sql.functions import monotonically_increasing_id

    # add unique identifier to column names
    #df = Z.toDF(*[f"{col}_{i}" for i, col in enumerate(Z.columns)])
    #df=df.toPandas()
    
    df=Z.iloc[:,1:]
    df=np.array(df)

    out=np.empty([len(list_of_ids),len(list_of_ids)])
    np.fill_diagonal(out, 1)
    i=1
    for i in range(0,len(list_of_ids)):
        for j in range(i+1,len(list_of_ids)):
            l=df[:,[i,j]]
            l=l[(l[:,0]**2<=4) & (l[:,1]**2<=4),]
            x=l[:,0]
            y=l[:,1]
            out[i,j]=np.corrcoef(x,y)[1,0]
            out[j,i]=out[i,j]

    spark.stop()

    return out


def GIP1_lin_comb_Z_based(Z, covm, eaf, N, gcor, h2):
    import numpy as np
    from scipy.stats import chi2
    import pandas as pd
    gcov=np.sqrt(np.outer(h2, h2))*gcor
    eigen_values=np.linalg.eig(gcov)[0].real
    a=np.linalg.eig(gcov)[1][:,np.argmax(eigen_values)].real
    if a[0]<0: a=-a
    vary = np.sum(covm * np.outer(a, a))
    a=a/np.sqrt(vary)
    vary = np.sum(covm * np.outer(a, a))
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
    SE = np.sqrt(1 / (N*(1+((Z**2)/N))))
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



def number_of_sig_snps(list_of_ids,max_ram_to_use):
    from pyspark.sql import SparkSession
    import pyspark.sql.functions as f
    import pyspark.sql.types as t
    from pyspark.sql import Window
    import scipy as sc
    from scipy import stats
    import numpy as np
    import pandas as pd

    #
    global spark
    spark = (
        SparkSession.builder
        .master('local[*]')
        .config('spark.driver.memory', max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )

    def create_empty_dataframe(rows, columns):
        df = pd.DataFrame(index=range(rows), columns=range(columns))
        return df


    out=create_empty_dataframe(rows=len(list_of_ids),columns=4)

    gw=list_of_ids[0]
    i=0
    for gw in list_of_ids:
        print(gw)
        gw_path='gs://genetics-portal-dev-sumstats/unfiltered/gwas/'+gw+'.parquet/'
        gwas=spark.read.parquet(gw_path)
        out.iloc[i,0]=gw
        gwas=(gwas
        .filter(f.col("eaf")>=0.001)
        .filter(f.col("eaf")<=0.999)
        )
        out.iloc[i,1]=gwas.count()
        gwas=(gwas
            .filter(f.col("pval")<=5e-8))
        out.iloc[i,2]=gwas.count()
        out.iloc[i,3]=out.iloc[i,2]/out.iloc[i,1]
        i=i+1
    
    spark.stop()
    return out

