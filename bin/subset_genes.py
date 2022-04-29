#!/usr/bin/env python3

import argparse
import pandas as pd
# import numpy as np
from math import nan
from itertools import compress
from collections import Counter

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gene_presence_absence",
        help = "CSV file containing presence/absence of genes per clone"
    )
    parser.add_argument(
        "--metadata_merged",
        help = "TSV file mapping clone names with sequence data sample names"
    )
    parser.add_argument(
        "--perc",
        help = "Minimum percentage of strains containing gene for inclusion (optional)",
        nargs = '?', default = None
    )
    parser.add_argument(
        "--rm_split", default = True,
        help = "Remove genes which are marked as split? Default = True",
        type = lambda x: (str(x).lower() == 'true')
    )
    parser.add_argument(
        "--ref_only", default = False,
        help = "Only keep genes that are present in the reference strain? Default = False",
        type = lambda x: (str(x).lower() == 'true')
    )
    parser.add_argument("--outf", help="File for results")
    args = parser.parse_args()
    find_core(**vars(args))

def find_core(gene_presence_absence, metadata_merged, perc, rm_split, ref_only, outf):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    if ref_only: ## only include genes present in the reference
        csv_data = csv_data[csv_data.iloc[:,3].notna()]
        csv_data = csv_data.reset_index(drop=True)
    metadata = pd.read_csv(metadata_merged, sep = "\t")
    colnames = csv_data.columns.values.tolist()
    gene_names = csv_data.iloc[:,0].tolist()
    ## subset metadata to strains present in the presence/absence CSV
    metadata = metadata[metadata['DNA_sample_id'].isin(colnames)]
    clone_names = metadata['DNA_sample_id'].tolist()
    clone_data = csv_data[csv_data.columns.intersection(clone_names)]
    clone_sub = clone_data[clone_names].copy(deep=True)
    colnames_clone = clone_data.columns.values.tolist()
    ## sort metadata to match the clone_data
    metadata.DNA_sample_id = metadata.DNA_sample_id.astype("category")
    metadata.DNA_sample_id.cat.set_categories(colnames_clone)
    metadata = metadata.sort_values(["DNA_sample_id"])
    if rm_split:
        split_count = (clone_sub.apply(lambda x: x.str.findall(';').str.len()))
        split_count = split_count.sum(axis=1).astype(int)
        # rf_counts = (clone_sub.apply(lambda x: x.str.findall('refound').str.len()))
        # rf_counts = rf_counts.sum(axis=1).astype(int)
        # total_counts = rf_counts.add(split_count)
        total_counts = split_count
        total_counts.name = "total_counts"
        rm_idx = list(total_counts.loc[total_counts>0].index.to_numpy())
        for index in sorted(rm_idx, reverse=True):
            # NB delete in reverse order to avoid throwing off the subsequent indexes.
            del gene_names[index]
        clone_sub.drop(clone_sub.index[rm_idx], inplace=True)
    else:
        clone_sub.replace(';', nan, regex=True, inplace=True)
        clone_sub.replace('refound', nan, regex=True, inplace=True)
    ## subset to genes present in at least `perc` strains
    na_counts = clone_sub.isnull().sum(axis=1)
    max_na = len(clone_sub.columns)-(round((len(clone_sub.columns)/100)*float(perc)))
    clone_sub['gene'] = gene_names
    clone_sub = clone_sub[na_counts <= max_na]
    clone_sub.to_csv(outf, index=False, sep='\t')



if __name__ == "__main__":
    parse()
