import numpy as np
import pandas as pd
from numpy.linalg import inv

def GWAS_linear_combination_v2(a, beta_a, se, var_y=None, covm=None, N=None):
    snp_row = np.arange(beta_a.shape[0])
    if var_y is None:
        var_y = np.ones(a.shape[0])
    
    for i in range(var_y.shape[0]):
        beta_a[:, i] = beta_a[:, i] / np.sqrt(var_y[i])
        se[:, i] = se[:, i] / np.sqrt(var_y[i])

    if covm is not None:
        covm = (1 / (np.sqrt(np.diag(covm))[:, np.newaxis] @ np.sqrt(np.diag(covm))[np.newaxis, :])) * covm

    var_y_a = np.sum(covm * (a[:, np.newaxis] @ a[np.newaxis, :]))
    b = beta_a @ a
    varb = np.array([np.sum(covm * (se[i, np.newaxis] @ se[i, np.newaxis].T) * (a[:, np.newaxis] @ a[np.newaxis, :])) for i in snp_row])
    sen = np.sqrt(varb)
    b = b / np.sqrt(var_y_a)
    sen = sen / np.sqrt(var_y_a)
    
    ind = np.where(np.abs(b / sen) < 2)[0]
    N_from_se = np.median(1 / (sen[ind]**2))
    
    out = pd.DataFrame({'b': b, 'se': sen, 'N': N_from_se})
    
    return out


def GWAS_linear_combination_Z_based(a, Z, covm, N, eaf):
    snp_row = np.arange(Z.shape[0])
    var_y = np.ones(a.shape[0])
    se_a = np.apply_along_axis(lambda x: np.sqrt(1 / (x[0]**2 + x[1])), axis=1, arr=np.column_stack((Z, N)))
    beta_a = np.apply_along_axis(lambda x: x[0] * se_a[x[2], :], axis=1, arr=np.column_stack((Z, N, snp_row)))
    
    out = GWAS_linear_combination_v2(a=a, beta_a=beta_a, se=se_a, var_y=None, covm=covm, N=N)
    
    varg = 2 * (1 - eaf) * eaf
    out['b'] = out['b'] / np.sqrt(varg)
    out['se'] = out['se'] / np.sqrt(varg)
    
    return out
