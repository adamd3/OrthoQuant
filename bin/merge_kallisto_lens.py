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
        "--quant_dir",
        help = "Directory containing Kallisto quantification files"
    )
    parser.add_argument(
        "--metadata_merged",
        help = "TSV file mapping clone names with sequence data sample names"
    )
    parser.add_argument(
        "--ST_file",
        help = "File containing specific STs to be included (optional)",
        nargs = '?', default = None
    )
    parser.add_argument("--outf", help="File for results")
    args = parser.parse_args()
    merge_lens(**vars(args))

def merge_lens(gene_presence_absence, quant_dir, metadata_merged, ST_file, outf):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    metadata = pd.read_csv(metadata_merged, sep = "\t")
    st_tab = pd.read_csv(ST_file, sep = "\t", header=None)
    colnames = csv_data.columns.values.tolist()
    gene_names = csv_data.iloc[:,0].tolist()
    metadata = metadata[metadata['sample_id'].isin(colnames)]
    ## subset to STs with specific list
    keep_ST = st_tab[0].tolist()
    keep_ST = [str(st) for st in keep_ST]
    meta_sub = metadata[metadata['majority_ST'].isin(keep_ST)]
    quant_dfs = []
    for index, row in meta_sub.iterrows():
        sample_id = row['sample_id']
        sample_name = row['sample_name']
        quant_file = os.path.join(quant_dir, 'abundance.tsv')
        quant_dat = pd.read_csv(quant_file, sep = "\t")
        quant_dat = quant_dat[["target_id", "eff_length"]]
        quant_dat = quant_dat.rename(columns={'target_id': sample_id})
        cg = csv_data[["Gene", sample_id]]
        quant_merged = pd.merge(cg, quant_dat, on=sample_id)
        quant_merged = quant_merged[["Gene", "eff_length"]]
        quant_merged = quant_merged.rename(columns={'eff_length': sample_name})
        quant_dfs.append(quant_merged)
    quant_dfs = [df.set_index('Gene') for df in quant_dfs]
    quant_merged = pd.concat(quant_dfs, axis=1)
    quant_merged.to_csv(outf, index=True, sep='\t')



if __name__ == "__main__":
    parse()
