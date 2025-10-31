import os
import time

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
    representatives_of_clusters={}
    for cl_no, prot_id in clusters.items():
        old_header = proteins[prot_id]["header"]
        new_header = old_header.replace(">", f">{cl_no};")
        new_clusters[prot_id] = {"header": new_header, "seq": proteins[prot_id]["seq"]}
        representatives_of_clusters[cl_no] = {"header": new_header, "seq": proteins[prot_id]["seq"]}
    return new_clusters, representatives_of_clusters


def save_clusters_representatives(file, new_proteins):
    with open(file, "w") as f:
        for prot, prot_info in new_proteins.items():
            f.write(prot_info["header"])
            f.write(prot_info["seq"])


def blast(new_proteins, folder):
    rids = []
    if not os.path.exists(folder):
        os.mkdir(folder)
    new_proteins_no = len(new_proteins)
    new_protein_data = list(new_proteins.values())
    while new_protein_data:
        new_protein_data_tmp = new_protein_data[:20]
        new_protein_data = new_protein_data[20:]
        protein = ''
        for new_protein_info in new_protein_data_tmp:
            protein += new_protein_info["header"] + new_protein_info["seq"]
        rids.append(blast_req(protein))
        print(len(rids), len(new_protein_data), rids)
    NEW_RIDS = []
    e = 0
    with open("tmp.csv", "w") as f:
        f.write(",".join([str(i) for i in rids]))
    while rids:
        URl = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rids[e]}"
        req = requests.get(URl)
        status = req.text.split("QBlastInfoBegin")[-1].split("QBlastInfoEnd")[0]
        if status.splitlines()[1].split("=")[-1] == "yes":
            url = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID={rids[e]}"
            req = requests.get(url)
            with open(os.path.join(folder, f"{rids[e]}.txt")) as f:
                f.write(req.text)
        else:
            NEW_RIDS.append(rids[e])
        if len(rids) < 10:
            time.sleep(10)
        e += 1
        if len(rids) == e:
            rids = NEW_RIDS
            print(len(rids))
            NEW_RIDS = []
            e = 0


def blast_req(protein):
    request = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr_cluster_seq",
        "QUERY": protein,
        "ALIGNMENTS": 1000
    }
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    req = requests.post(url, data=request)
    RID = req.text.split("QBlastInfoBegin")[1].split("QBlastInfoEnd")[0].splitlines()[1].split()[-1]
    return RID

def read_corr_file(result, file):
    with open(file) as f:
        for l in f:
            if l.strip():
                line=l.strip().split(";")
                metabolite=line[0].strip()
                pval=float(line[4])
                cluster=line[1].split("_")[1]
                if pval<0.05 and metabolite in {"Conc Cholic acid", "Conc Litocholic acid"}:
                    result.add(cluster)
    return result
                

def read_corr_files(files):
    result=set()
    for file in files:
        result=read_corr_file(result, file)
    return result

def select_representative_significant_sequences(proteins_with_corr,new_proteins_by_clusters):
    result={}
    for cl in list(proteins_with_corr):
        result[cl]=new_proteins_by_clusters[cl]
    return result

@click.command()
@click.option('--clusters_file', default="./", help='File with count of genes.')
@click.option('--representatives_of_clusters', default="./", help='File with count of genes.')
@click.option('--cluster_representatives_file', default="./", help='Out file with normalised genes statistics.')
@click.option('--cluster_representatives_with_corr_file', default="./", help='Out file with normalised genes statistics.')
@click.option('--out_folder', default="./", help='Out file with normalised genes statistics.')
@click.option('--corr_files', default="./", help='Out file with normalised genes statistics.')
def run_blast(clusters_file, representatives_of_clusters, cluster_representatives_file, out_folder, corr_files, cluster_representatives_with_corr_file):
    clusters = read_clusters_representatives(clusters_file)
    proteins = read_representative_fasta(representatives_of_clusters)
    new_proteins, new_proteins_by_clusters = update_proteins_names(proteins, clusters)
    save_clusters_representatives(cluster_representatives_file, new_proteins)
    proteins_with_corr = read_corr_files(corr_files.strip().split(','))
    representative_significant_sequences=select_representative_significant_sequences(proteins_with_corr,new_proteins_by_clusters)
    save_clusters_representatives(cluster_representatives_with_corr_file, representative_significant_sequences)
    #blast(representative_significant_sequences, out_folder)
    # makeblastdb -in ./uniprotkb_ec_3_5_1_24_2025_10_22.fasta -dbtype prot
    # makeblastdb -in ./uniprot_trembl.fasta -dbtype prot
    #  blastp -db ./uniprotkb_ec_3_5_1_24_2025_10_22.fasta -out ./3_5_1_24_proteins -query ./clusters_representatives_for_corr.fasta -num_threads 10 -outfmt 6 -evalue 0.01 -qcov_hsp_perc 90




if __name__ == "__main__":
    run_blast()
