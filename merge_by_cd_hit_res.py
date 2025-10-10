import os

import click
import numpy as np


def read_clusters_file(path):
    clusters = {}
    cl_name = ""
    cl_sequences = set()
    with open(path) as f:
        #         >Cluster 53
        # 0       4447aa, >5928034_9... *
        # 1       4242aa, >2370890_1... at 99.53%
        # 2       4447aa, >3628458_10... at 99.35%
        for l in f:
            if l.startswith(">"):
                clusters[cl_name] = cl_sequences
                cl_sequences = set()
                cl_name = l.strip().split()[-1]
            else:
                cl_sequences.add(l.split(">")[-1].split("...")[0])
        if cl_sequences:
            clusters[cl_name] = cl_sequences
    return clusters


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


def merge_by_clusters(clusters, stats_file):
    stats_file_tmp = {}
    for file, genes_data in stats_file.items():
        new_stats = {}
        stats_file_tmp[file] = {}
        for cluster, proteins_in_cl in clusters.items():
            proteins_in_cl = list(proteins_in_cl)
            if len(proteins_in_cl) == 1:
                if proteins_in_cl[0] in new_stats:
                    new_stats[proteins_in_cl[0]] = genes_data[proteins_in_cl[0]]
            else:
                for prot in proteins_in_cl:
                    if prot in stats_file[file]:
                        if f"cluster_{cluster}" not in new_stats:
                            new_stats[f"cluster_{cluster}"] = genes_data[prot]
                        else:
                            new_stats[f"cluster_{cluster}"] += genes_data[prot]
        stats_file_tmp[file] = new_stats
    return stats_file_tmp


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


@click.command()
@click.option('--count_genes_file', default="./", help='File with count of genes.')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
@click.option('--cd_hit_res', default="./", help='File with CD-HIT clusters.')
def merge(count_genes_file, out_file, cd_hit_res):
    clusters = read_clusters_file(cd_hit_res)
    genes_counts = read_stats_file(count_genes_file)
    merged_stats = merge_by_clusters(clusters, genes_counts)
    save_new_stats(merged_stats, out_file)


if __name__ == "__main__":
    merge()
