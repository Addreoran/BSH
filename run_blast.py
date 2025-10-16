import os

import click
import numpy as np
from scipy.stats import pearsonr
import requests


def read_clusters_representatives(file):
    cluster_no = None
    clusters = {}
    with open(file) as f:
        for l in f:
            if l.startswith(">"):
                cluster_no = l.strip().split()[-1]
            elif "*" in l:
                clusters[cluster_no] = l.strip().split(">")[-1].split("...")[0]
    return clusters


def read_representative_fasta(file):
    representatives = {}
    sequence = ''
    header = None
    protein = None
    with open(file) as f:
        for l in f:
            if l.startswith(">"):
                if header is not None:
                    representatives[protein] = {"header": header, "seq": sequence}
                    header = None
                    protein = None
                    sequence = ''
                header = l
                protein = l.split()[0][1:]
            else:
                sequence += l
        if header is not None:
            representatives[protein] = {"header": header, "seq": sequence}
    return representatives


def update_proteins_names(proteins, clusters):
    new_clusters = {}
    for cl_no, prot_id in clusters.items():
        old_header = proteins[prot_id]["header"]
        new_header = old_header.replace(">", f">{cl_no}")
        new_clusters[prot_id] = {"header": new_header, "seq": proteins[prot_id]["seq"]}
    return new_clusters


def save_clusters_representatives(file, new_proteins):
    with open(file, "w") as f:
        for prot, prot_info in new_proteins.items():
            f.write(prot_info["header"])
            f.write(prot_info["seq"])


def blast(file):
    request = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "",
        "QUERY": file.open().read()
    }
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    req = requests.post(url, data=request)
    print(req.content)


@click.command()
@click.option('--clusters_file', default="./", help='File with count of genes.')
@click.option('--representatives_of_clusters', default="./", help='File with count of genes.')
@click.option('--cluster_representatives_file', default="./", help='Out file with normalised genes statistics.')
@click.option('--out_file', default="./", help='Out file with normalised genes statistics.')
def run_blast(clusters_file, representatives_of_clusters, cluster_representatives_file, out_file):
    clusters = read_clusters_representatives(clusters_file)
    proteins = read_representative_fasta(representatives_of_clusters)
    new_proteins = update_proteins_names(proteins, clusters)
    save_clusters_representatives(cluster_representatives_file, new_proteins)
    blast(cluster_representatives_file)


if __name__ == "__main__":
    run_blast()
