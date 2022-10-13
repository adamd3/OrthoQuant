#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import os.path

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sample_file",
        help = "tsv file containing sample metadata"
    )
    parser.add_argument("outf", help="File for output")
    args = parser.parse_args()
    check_meta(**vars(args))

def check_meta(sample_file, outf):
    sample_dat = pd.read_csv(sample_file, sep = "\t")
    # remove empty rows
    sample_dat.replace("", float("NaN"), inplace=True)
    sample_dat.dropna(subset = ["dna_sample_id"], inplace=True)
    if(len(sample_dat.fastq2.value_counts()) > 0):
        sample_dat['paired'] = "1"
    else:
        sample_dat['paired'] = "0"
    sample_dat.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()
