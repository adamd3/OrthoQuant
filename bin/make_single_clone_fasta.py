#!/usr/bin/env python3

import argparse
import pandas as pd
import os.path
import textwrap
from Bio import SeqIO


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
    parser.add_argument(
        "strain_name",
        help = "Name of the strain for which FASTA file will be created"
    )
    args = parser.parse_args()
    make_fasta(**vars(args))

def make_fasta(multifasta_file, gene_presence_absence, strain_name):
    csv_data = pd.read_csv(gene_presence_absence, low_memory=False)
    colnames = csv_data.columns.values.tolist()
    # gene_seq = Fasta(multifasta_file)
    if strain_name in colnames:
        strain_genes = csv_data[strain_name].tolist()
        records = (r for r in SeqIO.parse(multifasta_file, "fasta") if r.id in strain_genes)
        # genes_present = [gene for gene in strain_genes if gene in gene_seq.keys()]
        outf2 = os.path.normpath(strain_name+'.fna')
        SeqIO.write(records, outf2, "fasta")
        # f = open(outf2, 'a')
        # for gene in genes_present:
        #     seq = '\n'.join(textwrap.wrap(str(gene_seq[gene]), 50))
        #     f.write(">"+gene+"\n"+seq+"\n")
        # f.close()

if __name__ == "__main__":
    parse()
