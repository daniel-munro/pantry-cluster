"""Assemble data into an RNA phenotype BED file"""

import argparse
from pathlib import Path
from gtfparse import read_gtf
import numpy as np
import pandas as pd

def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations from GTF file
    
    Returns TSS as the first four columns of the BED format, meaning the
    coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'gene', :]
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    anno = anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]
    # Rename columns for tensorQTL:
    anno.columns = ['#chr', 'start', 'end', 'gene_id']
    return anno

def transcript_to_gene_map(ref_anno: Path) -> pd.DataFrame:
    """Load transcript IDs and corresponding gene IDs from GTF file"""
    anno = read_gtf(ref_anno)
    anno = anno.loc[anno['feature'] == 'transcript', :]
    return anno[['gene_id', 'transcript_id']]

def load_kallisto(sample_ids: list, kallisto_dir: Path, units: str) -> pd.DataFrame:
    """Assemble kallisto est_counts or tpm outputs into a table"""
    counts = []
    for i, sample in enumerate(sample_ids):
        fname = kallisto_dir / sample / 'abundance.tsv'
        d = pd.read_csv(fname, sep='\t', index_col='target_id')
        d = d[units] # est_counts or tpm
        d.name = sample
        # Store in lists and concat all at once to avoid 'PerformanceWarning: DataFrame is highly fragmented' warning
        counts.append(d)
    return pd.concat(counts, axis=1)

parser = argparse.ArgumentParser(description='Assemble kallisto est_counts into a gene-level BED file')
parser.add_argument('--input-dir', type=Path, required=True, help='Directory containing input files, for phenotype groups with per-sample input files')
parser.add_argument('--samples', type=Path, required=False, help='File listing sample IDs, for phenotype groups with per-sample input files')
parser.add_argument('--ref-anno', type=Path, required=True, help='Reference annotation file')
parser.add_argument('--output', type=Path, required=True, help='Output file ("*.bed")')
args = parser.parse_args()

if args.samples is not None:
    samples = pd.read_csv(args.samples, sep='\t', header=None)[0].tolist()

df_iso = load_kallisto(samples, args.input_dir, 'est_counts')
df_iso.index = df_iso.index.rename('transcript_id')
gene_map = transcript_to_gene_map(args.ref_anno).set_index('transcript_id')
df_iso = df_iso.join(gene_map, how='inner')
assert df_iso.shape[0] > 0, 'No matching isoforms in kallisto output and reference annotation'
# Sum to gene-level expression:
df_gene = df_iso.reset_index().drop('transcript_id', axis=1).groupby('gene_id', group_keys=True).sum()

anno = load_tss(args.ref_anno)
df_gene = anno.merge(df_gene.reset_index(), on='gene_id', how='inner')
df_gene = df_gene.rename(columns={'gene_id': 'phenotype_id'})
df_gene = df_gene[['#chr', 'start', 'end', 'phenotype_id'] + samples]
df_gene.to_csv(args.output, sep='\t', index=False, float_format='%g')
