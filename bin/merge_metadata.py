import argparse
import pandas as pd
import os.path

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "metadata_file", help = "contains clone ID and SRR ID for WGS data"
    )
    parser.add_argument(
        "data_dir", help = "Directory containing downloaded RNA-seq data"
    )
    parser.add_argument("outf", help="File for results")
    args = parser.parse_args()
    merge_meta(**vars(args))

def merge_meta(metadata_file, data_dir, outf):
    id_mappings = os.path.join(data_dir, 'samplesheet/id_mappings.csv')
    id_map_dat = pd.read_csv(id_mappings).iloc[:,[0,7]]
    clone_metadat = pd.read_csv(metadata_file, sep = "\t").iloc[:,:16]
    id_map_dat = id_map_dat.rename(
        columns={'sample': 'RNA_sample_id', 'sample_title': 'sample_name'}
    )
    merged_meta = clone_metadat.merge(id_map_dat, on='sample_name')
    ## move RNA column to front of df
    merged_meta = merged_meta[
        ['RNA_sample_id'] +
        [ col for col in merged_meta.columns if col != 'RNA_sample_id' ]
    ]
    merged_meta['fastq'] = [
        os.path.join(data_dir, 'fastq', s+'_T1.fastq.gz') for s in merged_meta['RNA_sample_id'].tolist()]
    merged_meta.to_csv(outf, index=False, sep='\t')

if __name__ == "__main__":
    parse()
