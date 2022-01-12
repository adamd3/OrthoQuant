#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import textwrap
from pyfaidx import Fasta

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "multifasta_file",
        help = "Multifasta file containing gene sequences"
    )
    parser.add_argument(
        "gene_presence_absence",
        help = "CSV file containing presence/absence of genes per clone"
    )
    parser.add_argument("outdir", help="Directory for results")
    args = parser.parse_args()
    make_fasta(**vars(args))

def make_fasta(multifasta_file, gene_presence_absence, outdir):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    colnames = csv_data.columns.values.tolist()
    gene_seq = Fasta(multifasta_file)
    for idx in range(4, len(csv_data.columns)):
        clone_name = colnames[idx]
        genes = csv_data.iloc[:,idx].tolist()
        genes_present = [gene for gene in genes if gene in gene_seq.keys()]
        if (len(genes_present)>0):
            outf = os.path.join(outdir, clone_name+'.fna')
            f = open(outf, 'a')
            for gene in genes_present:
                seq = '\n'.join(textwrap.wrap(str(gene_seq[gene]), 50))
                f.write(">"+gene+"\n"+seq+"\n")
            f.close()

if __name__ == "__main__":
    parse()
