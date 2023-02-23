import numpy as np
from scipy.stats import chi2
import pandas as pd

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

