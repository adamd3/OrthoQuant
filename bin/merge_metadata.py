#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "metadata_file", help = "contains clone ID and SRR ID for WGS data"
    )
    parser.add_argument(
        "sample_ID_file", help = "contains DNA > RNA sample ID mappings"
    )
    parser.add_argument(
        "data_dir", help = "Directory containing downloaded RNA-seq data"
    )
    parser.add_argument("outf", help="File for results")
    args = parser.parse_args()
    merge_meta(**vars(args))

def merge_meta(metadata_file, sample_ID_file, data_dir, outf):
    id_map_dat = pd.read_csv(sample_ID_file).iloc[:,[0,7]]
    clone_metadat = pd.read_csv(metadata_file, sep = "\t").iloc[:,:16]
    id_map_dat = id_map_dat.rename(
        columns={
        'sample': 'rna_sample_id',
        'sample_title': 'sample_name'
        }
    )
    merged_meta = clone_metadat.merge(id_map_dat, on='sample_name')
    # merged_meta = merged_meta.rename(columns={'sample_id': 'dna_sample_id'})
    ## move RNA column to front of df
    merged_meta = merged_meta[['rna_sample_id'] +
        [ col for col in merged_meta.columns if col != 'rna_sample_id' ]
    ]
    # merged_meta['fastq'] = [os.path.join(data_dir, 'fastq', s +
    #     '_T1.fastq.gz') for s in merged_meta['rna_sample_id'].tolist()]
    # merged_meta['fasta'] = [os.path.join(data_dir, 'extract_fasta', s +
    #     '.fna') for s in merged_meta['dna_sample_id'].tolist()]
    ## remove any row where the FastQ or FASTA data is missing
    rm_idx = []
    for index, row in merged_meta.iterrows():
        if not os.path.exists(row['fastq']):
            rm_idx.append(index)
        if not os.path.exists(row['fasta']):
            rm_idx.append(index)
    merged_meta = merged_meta.drop(merged_meta.index[rm_idx])
    merged_meta.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()
