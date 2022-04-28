#!/usr/sbin/python3

import argparse
# import pandas as pd
import os.path
import csv
# import textwrap
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
    parser.add_argument(
        "strain_name",
        help = "Name of the strain for which FASTA file will be created"
    )
    args = parser.parse_args()
    make_fasta(**vars(args))

def make_fasta(multifasta_file, gene_presence_absence, strain_name):
    input_file = csv.DictReader(open(gene_presence_absence))
    strain_genes = []
    for row in input_file:
        gene = (row[strain_name])
        strain_genes.append(gene)
    strain_genes = list(filter(None, strain_genes))
    gene_seq = Fasta(multifasta_file)
    genes_present = [gene for gene in strain_genes if gene in gene_seq.keys()]
    outf = os.path.normpath(strain_name+'_ss.fna')
    f = open(outf, 'a')
    for gene in genes_present:
        seq = str(gene_seq[gene])
        f.write(">"+gene+"\n"+seq+"\n")
    f.close()

if __name__ == "__main__":
    parse()
