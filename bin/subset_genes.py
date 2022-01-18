#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
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
        "--ST_file",
        help = "File containing specific STs to be included",
        nargs = '?', default = None
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
    parser.add_argument("--outf", help="File for results")
    args = parser.parse_args()
    find_core(**vars(args))

def find_core(gene_presence_absence, metadata_merged, ST_file, perc, rm_split, outf):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    metadata = pd.read_csv(metadata_merged, sep = "\t")
    st_tab = pd.read_csv(ST_file, sep = "\t", header=None)
    colnames = csv_data.columns.values.tolist()
    gene_names = csv_data.iloc[:,0].tolist()
    ## subset metadata to strains present in the presence/absence CSV
    metadata = metadata[metadata['sample_id'].isin(colnames)]
    clone_names = metadata['sample_id'].tolist()
    clone_data = csv_data[csv_data.columns.intersection(clone_names)]
    colnames_clone = clone_data.columns.values.tolist()
    ## sort metadata to match the clone_data
    metadata.sample_id = metadata.sample_id.astype("category")
    metadata.sample_id.cat.set_categories(colnames_clone)
    metadata = metadata.sort_values(["sample_id"])
    ## Subset sequence types
    keep_ST = st_tab[0].tolist()
    keep_ST = [str(st) for st in keep_ST]
    meta_sub = metadata[metadata['majority_ST'].isin(keep_ST)]
    clone_sub = clone_data[meta_sub['sample_id'].tolist()].copy(deep=True)
    if rm_split:
        split_count = (clone_sub.apply(lambda x: x.str.findall(';').str.len()))
        split_count = split_count.sum(axis=1).astype(int)
        total_counts = split_count
        total_counts.name = "total_counts"
        rm_idx = list(total_counts.loc[total_counts>0].index.to_numpy())
        for index in sorted(rm_idx, reverse=True):
            # NB delete in reverse order to avoid throwing off the subsequent indices.
            del gene_names[index]
        clone_sub.drop(clone_sub.index[rm_idx], inplace=True)
    else:
        clone_sub.replace(';', np.NaN, regex=True, inplace=True) ## split genes
        clone_sub.replace('refound', np.NaN, regex=True, inplace=True) ## refound genes
    if perc is not None:
        ## subset to genes present in at least `perc` strains
        na_counts = clone_sub.isnull().sum(axis=1)
        max_na = len(clone_sub.columns)-(round((len(clone_sub.columns)/100)*float(perc)))
        clone_sub['gene'] = gene_names
        clone_sub = clone_sub[na_counts < max_na]
        clone_sub.to_csv(outf, index=False, sep='\t')
    else:
        ## get the core genome for each included ST
        clone_data_t = clone_sub.transpose()
        clone_data_t = clone_data_t.notnull().astype('int')
        clone_data_t.columns = gene_names
        clone_data_t['majority_ST'] = meta_sub['majority_ST'].tolist()
        group_dfs = dict(tuple(clone_data_t.groupby('majority_ST')))
        group_cores = {}
        for group in group_dfs.keys():
            df = group_dfs[group]
            df = df.drop(['majority_ST'], axis=1)
            nrow = len(df.index)
            colcounts = df.sum(axis=0)  ## count presence per gene
            keepcols = (colcounts == nrow)
            df = df[keepcols.index[keepcols]]
            group_cores[group] = df.columns.values.tolist()
        core_merge = pd.DataFrame.from_dict(group_cores, orient='index').T
        core_merge.to_csv(outf, index=False, sep='\t')



if __name__ == "__main__":
    parse()
