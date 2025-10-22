import os

import click
import numpy as np

def read_corr(corr_file):
    pass

def compare_cntrl_corr(corr_file, cntrl_corr_info):
    pass

def get_database_sequences_info(fasta_database):
    pass

def read_blast_result(blast_table):
    pass

def fix_corr(corr_info, cntrl_corr_info):
    pass

def save_corr(fixed_corr, out_file):
    pass


@click.command()
@click.option('--corr_file', default="./", help='Folder with BLAST files.')
@click.option('--corr_ctrl_file', default="./", help='Folder with BLAST files.')
@click.option('--blast_table', default="./", help='')
@click.option('--fasta_database', default="./", help='')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
def main(corr_file, corr_ctrl_file, blast_table, fasta_database, out_file):
    corr_info = read_corr(corr_file)
    cntrl_corr_info = read_corr(corr_cntrl_file)
    blast_result = read_blast_result(blast_table)
    database_fasta_info = get_database_sequences_info(fasta_database)
    fixed_corr = fix_corr(corr_info, cntrl_corr_info)
    save_corr(fixed_corr, out_file)


if __name__=='__main__':
    main()
