import os

import click
import numpy as np
from scipy.stats import pearsonr


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
                        stats[file][gene] = float(genes_no[e])
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


def read_metabolites(count_metabolites_files):
    res = {}
    with open(count_metabolites_files) as f:
        header = f.readline()
        header = header.strip().split(";")[2:]
        for l in f:
            if l.strip():
                line = l.strip().split(";")
                probe_name = line[0]
            res[probe_name] = {}
            res[probe_name]["metabolites"] = {i: float(line[e + 2]) for e, i in enumerate(header)}
    return res


def select_files_of_metabolites(genes_counts, metabolites):
    for file, gene_data in genes_counts.items():
        # mgshot_S4358Nr27.1.fastq.gz.shotgun.tsv
        file_name = file.split("_", 1)[1].split["."][0]
        file_id = file.split("_", 1)[1].split["."][1]
        if file.split("_", 1)[1].split["."][0] in metabolites:
            metabolites[file_name][file_id] = gene_data
    return metabolites


def count_pearson_corr(metabolites):
    all_metabolites = set()
    all_genes = set()
    all_files = set()

    for files, metabolites_data in metabolites.items():
        all_metabolites = all_metabolites.union(set(metabolites_data["metabolites"].keys()))
        all_genes = all_genes.union(set(metabolites_data["1"].keys()))
        all_genes = all_genes.union(set(metabolites_data["2"].keys()))
        all_files.add(files)
    all_files = list(all_files)
    tests = {}  # (metabo, gene, no):{metabo_values:[], gene_values:[], pval:[], corr_value:[]}
    metabolites_sets = {m: [] for m in list(all_metabolites)}
    genes_1 = {m: [] for m in list(all_genes)}
    genes_2 = {m: [] for m in list(all_genes)}
    for file in all_files:
        for metabolite, metabolites_data in metabolites[file]["metabolites"].items():
            metabolites_sets[metabolite].append(metabolites_data)
        for gene_name, gene_no in metabolites[file]["1"].items():
            genes_1[gene_name].append(gene_no)
        for gene_name, gene_no in metabolites[file]["2"].items():
            genes_2[gene_name].append(gene_no)
    for metabolite, metabolites_sets in metabolites_sets.items():
        for gene_name, gene_sets in genes_1.items():
            correlation = pearsonr(metabolites_sets, gene_sets)
            tests[(metabolite, gene_name, 1)] = {"metabo_values": metabolites_sets, "gene_values": gene_sets,
                                                 "pval": correlation.pvalue, "corr_value": correlation.statistics}
        for gene_name, gene_sets in genes_2.items():
            correlation = pearsonr(metabolites_sets, gene_sets)
            tests[(metabolite, gene_name, 2)] = {"metabo_values": metabolites_sets, "gene_values": gene_sets,
                                                 "pval": correlation.pvalue, "corr_value": correlation.statistics}
    return tests


def save_corr_stats(corr, out_file):
    with open(out_file, "w") as f:
        for datasets_info, datasets_data in corr.items():
            f.write(';'.join(list(datasets_info)))
            f.write(";")
            f.write(str(datasets_data['corr_value']))
            f.write(";")
            f.write(str(datasets_data['pval']))
            f.write(";")
            f.write(','.join(datasets_data['metabo_values']))
            f.write(";")
            f.write(','.join(datasets_data['gene_values']))
            f.write("\n")


@click.command()
@click.option('--count_genes_file', default="./", help='File with count of genes.')
@click.option('--count_metabolites_files', default="./", help='File with count of genes.')
@click.option('--out_file', default="./", help='Out file with normalised genes statistics.')
def count_corr(count_genes_file, out_file, count_metabolites_files):
    genes_counts = read_stats_file(count_genes_file)
    metabolites = read_metabolites(count_metabolites_files)
    genes_of_metabolites = select_files_of_metabolites(genes_counts, metabolites)
    corr = count_pearson_corr(genes_of_metabolites)
    save_corr_stats(corr, out_file)


if __name__ == "__main__":
    count_corr()
