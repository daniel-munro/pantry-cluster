# Run in python script to use random_tiebreak option
# Based on examples at https://github.com/broadinstitute/tensorqtl

import argparse
import pandas as pd
from rpy2.robjects.packages import importr
import tensorqtl
from tensorqtl import genotypeio, cis, trans

parser = argparse.ArgumentParser(description="Run tensorQTL from command line but with extra options")
parser.add_argument("geno_prefix")
parser.add_argument("expression")
parser.add_argument("output", help="Output file path, or, for cis_nominal mode, output files prefix (with no directories)")
parser.add_argument("--covariates", help="Covariates file")
parser.add_argument("--cis_output", help="For cis_independent mode, output of tensorQTL in cis mode")
parser.add_argument("--groups", help="File with phenotype groups if applicable")
parser.add_argument("--window", type=int, default=1000000, help="cis-window size, default 1000000")
parser.add_argument("--output_dir", help="For cis_nominal mode, output directory")
parser.add_argument("--mode", required=True, help="Run mode: currently cis, cis_independent, cis_nominal, or genome_wide. genome_wide runs `map_trans` without removing cis-window variants.")
args = parser.parse_args()

# Check for rpy2 and qvalue packages first, since otherwise tensorQTL will do
# all the eQTL mapping and then fail.
_ = importr("qvalue")

if args.output_dir is not None:
    assert args.mode == 'cis_nominal'

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
chroms = variant_df['chrom'].unique()
keep = pheno_pos['chr'].isin(chroms)
pheno = pheno.loc[keep]
pheno_pos = pheno_pos.loc[keep]

if args.mode == "cis":
    d = cis.map_cis(genotype_df, variant_df, pheno, pheno_pos, covar, group_s=groups, window=args.window, random_tiebreak=True)
    tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    d.to_csv(args.output, sep="\t", float_format="%.6g")
elif args.mode == "cis_independent":
    cis_df = pd.read_csv(args.cis_output, sep="\t", index_col=0)
    d = cis.map_independent(genotype_df, variant_df, cis_df, pheno, pheno_pos, covar, group_s=groups, window=args.window, random_tiebreak=True)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
elif args.mode == "genome_wide":
    d = trans.map_trans(genotype_df, pheno, covariates_df=covar, maf_threshold=0) # cis has default min MAF 0, trans has default 0.05
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
elif args.mode == "cis_nominal":
    assert args.output_dir is not None
    assert args.groups is None
    cis.map_nominal(genotype_df, variant_df, pheno, pheno_pos, args.output, covariates_df=covar,
                    group_s=None, window=args.window, output_dir=args.output_dir)
else:
    print("Mode not recognized.")
