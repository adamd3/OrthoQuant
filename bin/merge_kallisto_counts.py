#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
from collections import Counter


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gene_presence_absence",
        help = "CSV file containing binary presence/absence of genes per clone"
    )
    parser.add_argument(
        "--metadata_merged",
        help = "TSV file mapping clone names with sequence data sample names"
    )
    parser.add_argument("--outf", help="File for results")
    args = parser.parse_args()
    merge_counts(**vars(args))

def merge_counts(gene_presence_absence, metadata_merged, outf):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    metadata = pd.read_csv(metadata_merged, sep = "\t")
    colnames = csv_data.columns.values.tolist()
    gene_names = csv_data.iloc[:,0].tolist()
    metadata = metadata[metadata['DNA_sample_id'].isin(colnames)]
    quant_dfs = []
    for index, row in metadata.iterrows():
        DNA_sample_id = row['DNA_sample_id']
        sample_name = row['sample_name']
        quant_file = os.path.join('kallisto_'+DNA_sample_id, 'abundance.tsv')
        quant_dat = pd.read_csv(quant_file, sep = "\t")
        quant_dat = quant_dat[["target_id", "est_counts"]]
        quant_dat = quant_dat.rename(columns={'target_id': DNA_sample_id})
        cg = csv_data[["Gene", DNA_sample_id]]
        quant_merged = pd.merge(cg, quant_dat, on=DNA_sample_id)
        quant_merged = quant_merged[["Gene", "est_counts"]]
        quant_merged = quant_merged.rename(columns={'est_counts': sample_name})
        quant_dfs.append(quant_merged)
    quant_dfs = [df.set_index('Gene') for df in quant_dfs]
    quant_merged = pd.concat(quant_dfs, axis=1)
    quant_merged.to_csv(outf, index=True, sep='\t')


if __name__ == "__main__":
    parse()
