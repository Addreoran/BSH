import os

import click
import numpy as np


def read_stats_file(count_genes_file):
    file_names = list()
    stats = {}
    with open(count_genes_file) as f:
        for l in f:
            if l.strip():
                if "genes" in l:
                    file_names = l.strip().split(";")[1:]
                    stats = {i: {} for i in file_names}
                else:
                    line = l.strip().split(";")
                    gene = line[0]
                    genes_no = line[1:]
                    for e, file in enumerate(file_names):
                        stats[file][gene] = int(genes_no[e])
    return stats


def save_new_stats(merged_stats, out_file):
    file_list = list(merged_stats.keys())
    genes = set()
    for file, genes_data in merged_stats.items():
        genes = genes.union(set(genes_data.keys()))
    with open(out_file, "w") as f:
        f.write("genes;")
        f.write(";".join(file_list))
        f.write("\n")
        for gene in list(genes):
            f.write(f"{str(gene)};")
            for file_name in file_list:
                f.write(f"{str(merged_stats[file_name].get(gene, 0))};")
            f.write("\n")


def count_normalisation(gene_counts):
    reads_no = {}
    for file, genes_no in gene_counts.items():
        if file not in reads_no:
            reads_no[file] = 0
        for gene, no in genes_no.items():
            reads_no[file] += no
    res = {file: {} for file in gene_counts.keys()}
    for file, genes_no in gene_counts.items():
        for gene, no in genes_no.items():
            res[file][gene] = no / reads_no[file]
    return res


@click.command()
@click.option('--count_genes_file', default="./", help='File with count of genes.')
@click.option('--out_file', default="./", help='Out file with normalised genes statistics.')
def normalize(count_genes_file, out_file):
    genes_counts = read_stats_file(count_genes_file)
    normalised_counts = count_normalisation(genes_counts)
    save_new_stats(normalised_counts, out_file)


if __name__ == "__main__":
    normalize()
