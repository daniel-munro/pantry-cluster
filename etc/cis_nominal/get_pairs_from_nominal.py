import argparse
from fastparquet import ParquetFile
import pandas as pd

CHROMS = [f'chr{x}' for x in range(1, 23)]

parser = argparse.ArgumentParser(description='Extract all pairs for a gene from full tensorQTL cis output.')
parser.add_argument('--prefix', help='Path prefix (including directories) for nominal parquet files ("{prefix}.{chrom}.parquet").')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--gene', help='Gene ID of phenotypes to extract. Assumed to be the prefix of the `phenotype_id`. All associations for these phenotypes will be extracted.')
group.add_argument('--cis-qtl', help='Output of tensorqtl cis mode, including columns `phenotype_id` and `variant_id`. Matching pairs will be extracted.')
parser.add_argument('--chrom', help='Load and extract from only one chromosome, either for efficiency (`--gene`) or filtering (`--cis-qtl`).')
parser.add_argument('--simplify-splice-ids', action='store_true', help='To match leafcutter-style splice junction phenotype IDs from different datasets, remove the last part, e.g. `:clu_35321_-` and any `chr` in chromosome names.')
parser.add_argument('--output', help='output file (TSV, can be *.gz)')
args = parser.parse_args()

if args.cis_qtl is not None:
    pairs = pd.read_csv(args.cis_qtl, sep="\t")
    if args.simplify_splice_ids:
        pairs['phenotype_id'] = pairs['phenotype_id'].str.replace(r':clu_.+$', '', regex=True)
        pairs['phenotype_id'] = pairs['phenotype_id'].str.replace('chr', '')
    # Save phenotype-variant pairs to retrieve per chromosome for efficient lookup
    pairs['chrom_id'] = pairs['variant_id'].str.split('_').str[0]
    pairs = pairs[['phenotype_id', 'variant_id', 'chrom_id']]
    pairs = pairs.apply(lambda row: (row['phenotype_id'], row['variant_id']), axis=1).groupby(pairs['chrom_id']).apply(set).to_dict()

matches = []
for chrom in CHROMS:
    if args.chrom is not None and chrom != args.chrom:
        assert args.chrom in CHROMS
        continue
    print(f"Processing chromosome {chrom}", flush=True)
    df = ParquetFile(f"{args.prefix}.{chrom}.parquet").to_pandas()
    if args.gene is not None:
        df = df.loc[df['phenotype_id'].str.startswith(args.gene), :]
    elif args.simplify_splice_ids is None:
        df = df.loc[df.apply(lambda row: (row['phenotype_id'], row['variant_id']) in pairs[chrom], axis=1), :]
    else:
        df['phenotype_id2'] = df['phenotype_id'].str.replace(r':clu_.+$', '', regex=True)
        df['phenotype_id2'] = df['phenotype_id2'].str.replace('chr', '')
        df = df.loc[df.apply(lambda row: (row['phenotype_id2'], row['variant_id']) in pairs[chrom], axis=1), :]
        df = df.drop(columns=['phenotype_id2'])
    matches.append(df)

matches = pd.concat(matches, axis=0)
matches.to_csv(args.output, sep="\t", index=False, float_format="%g")
print(f"Saved {len(matches)} pairs to {args.output}", flush=True)
