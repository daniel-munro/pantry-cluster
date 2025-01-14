# Run in python script to use random_tiebreak option
# Based on examples at https://github.com/broadinstitute/tensorqtl

import argparse
import os
import numpy as np
import pandas as pd
from rpy2.robjects.packages import importr
import tensorqtl
from tensorqtl import genotypeio, cis
import torch

parser = argparse.ArgumentParser(description="Run tensorQTL from command line but with extra options")
parser.add_argument("geno_prefix")
parser.add_argument("expression")
parser.add_argument("output", help="Output file")
parser.add_argument("--covariates", help="Covariates file")
parser.add_argument("--cis_output", help="For cis_independent mode, output of tensorQTL in cis mode")
parser.add_argument("--groups", required=True, help="File with phenotype groups if applicable")
parser.add_argument("--mode", required=True, help="Run mode: currently either cis or cis_independent")
parser.add_argument("--chrom", required=True, help="If supplied, only phenotypes for this chromosome will be run, and calculate_qvalues will not be run. Afterwards, concatenate all chromosome results and run that.")
parser.add_argument("--batch", required=True, type=int, help="If supplied, after filtering phenotypes to one chromosome, phenotypes (or groups if provided) are further split into 10 consecutive batches, named 0-9, and only those will be run.")
parser.add_argument("--subbatch", required=True, type=int, help="Run only subbatch N of 10 within the batch.")
args = parser.parse_args()

# Check for rpy2 and qvalue packages first, since otherwise tensorQTL will do
# all the eQTL mapping and then fail.
_ = importr("qvalue")

if torch.cuda.is_available():
    print(f'  * using GPU ({torch.cuda.get_device_name(torch.cuda.current_device())})')
else:
    print('  * WARNING: using CPU!')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.expression)
if args.covariates is not None:
    covar = pd.read_csv(args.covariates, sep="\t", index_col=0).T
else:
    covar = None
pr = genotypeio.PlinkReader(args.geno_prefix)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
if args.groups is not None:
    groups = pd.read_csv(args.groups, sep="\t", index_col=0, header=None).squeeze('columns')
else:
    groups = None

# Remove phenotypes on chromosomes not present in genotypes to avoid error
if args.chrom is None:
    chroms = variant_df['chrom'].unique()
else:
    chroms = [args.chrom]
keep = pheno_pos['chr'].isin(chroms)
pheno = pheno.loc[keep]
pheno_pos = pheno_pos.loc[keep]

# Filter groups to those with phenotypes on this chromosome
groups = groups.loc[groups.index.isin(pheno.index)]

# Assign each group to one of 10 batches of consecutive groups, in current order
group_ids = groups.unique()
batch_size = (len(group_ids) + 10 - 1) // 10
group_ids = [g for i, g in enumerate(group_ids) if i // batch_size == args.batch]

# Filter to subbatch
subbatch_size = (len(group_ids) + 10 - 1) // 10
group_ids = [g for i, g in enumerate(group_ids) if i // subbatch_size == args.subbatch]

# Get phenotypes in these groups
groups = groups.loc[groups.isin(group_ids)]
keep = pheno.index.isin(groups.index)
pheno = pheno.loc[keep]
pheno_pos = pheno_pos.loc[keep]

if args.mode == "cis":
    d = cis.map_cis(genotype_df, variant_df, pheno, pheno_pos, covar, group_s=groups, random_tiebreak=True)
    if args.chrom is None:
        tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    d.to_csv(args.output, sep="\t", float_format="%.6g")
elif args.mode == "cis_independent":
    cis_df = pd.read_csv(args.cis_output, sep="\t", index_col=0)
    cis_df = cis_df.loc[cis_df.index.isin(pheno.index)]
    d = cis.map_independent(genotype_df, variant_df, cis_df, pheno, pheno_pos, covar, group_s=groups, random_tiebreak=True)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
else:
    print("Mode not recognized.")
