"""Load expression read count BED file, log transform and normalize, regress out covariates, and save to CSV for aFC-n"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

def log_and_normalize(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Log-transform and normalize gene expression matrix
    
    Rows are genes. Adapted from aFCn code:
    https://github.com/PejLab/aFCn/blob/main/src/calc.pyx#L96
    """
    counts = counts_df.to_numpy()
    original_counts = counts.copy()
    original_counts = original_counts + 1
    counts = counts + 1
    # counts = counts[np.alltrue(counts > 100, axis=1)]
    logcounts = np.log2(counts)
    log_original_counts = np.log2(original_counts)
    # median of rows
    loggeommeans = np.median(logcounts, axis=1).reshape(len(logcounts), 1)
    loggeommeans_expand = np.repeat(loggeommeans, logcounts.shape[1], axis=1)
    # median of columns
    sf = np.median(logcounts - loggeommeans_expand, axis=0).reshape(1, logcounts.shape[1])
    sf_expand = np.repeat(sf,log_original_counts.shape[0], axis=0)
    normalized_matrix = log_original_counts - sf_expand
    # Convert back to pandas DataFrame
    df = pd.DataFrame(normalized_matrix, index=counts_df.index, columns=counts_df.columns)
    return df

def regress_out_covars(y: pd.Series, covars: np.array) -> pd.Series:
    """Identify strongly correlated covariates and regress them out"""
    # if y.name[-2:] == '00':
    #     print(y.name, flush=True)
    # Identify covariates with p < 0.01
    cov = sm.add_constant(covars)
    model1 = sm.OLS(y, cov).fit()
    assert model1.pvalues.index[0] == 'const'
    covars_sig = covars[:, model1.pvalues[1:] < 0.01] # exclude intercept
    if covars_sig.shape[1] == 0:
        return y
    # Regress out significant covariates
    model2 = LinearRegression().fit(covars_sig, y)
    resid = y - model2.predict(covars_sig)
    return pd.Series(resid, index=y.index)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-b', '--bed', type=Path, help='Expression BED file, non-log-transformed read counts')
parser.add_argument('-c', '--covar', type=Path, help='Phenotypes CSV file')
parser.add_argument('-o', '--output', type=Path, help='Output CSV file')
args = parser.parse_args()

# Load expression BED file
df = pd.read_csv(args.bed, sep='\t', index_col=3, dtype={'#chr': str, 'start': int, 'end': int})
df = df.drop(['#chr', 'start', 'end'], axis=1)

covar = pd.read_csv(args.covar, sep='\t', index_col=0)
assert df.columns.equals(covar.columns), 'Phenotype names do not match between expression and covariates'
covar = covar.to_numpy().T

# # Test
# df = df.iloc[:100, :]

df = log_and_normalize(df)
# Regress out covariates and save
df = df.apply(regress_out_covars, axis=1, args=(covar,)) # axis=1 means by row, i.e. index of Series is column names
df.index.name = 'Name'
df.to_csv(args.output)
