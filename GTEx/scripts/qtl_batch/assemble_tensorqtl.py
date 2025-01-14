import argparse
import os
import pandas as pd
from rpy2.robjects.packages import importr
import tensorqtl

parser = argparse.ArgumentParser(description="Assemble tensorqtl results from all chromosomes and do qvalue etc.")
parser.add_argument("input_prefix", help="Path prefix to inputs, e.g. the 'path/modality.cis_qtl' part of 'path/modality.cis_qtl.chr1.txt.gz'")
parser.add_argument("output", help="Output file")
parser.add_argument("--mode", required=True, help="Run mode: currently either cis or cis_independent")
args = parser.parse_args()

if args.mode == "cis":
    d = [pd.read_csv(f"{args.input_prefix}.chr{i}.txt.gz", sep="\t", index_col=0) for i in range(1, 23)]
    d = pd.concat(d)
    tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    d.to_csv(args.output, sep="\t", float_format="%.6g")
elif args.mode == "cis_independent":
    d = [pd.read_csv(f"{args.input_prefix}.chr{i}.txt.gz", sep="\t") for i in range(1, 23)]
    d = pd.concat(d)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
else:
    print("Mode not recognized.")
