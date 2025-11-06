import os

import click
import numpy as np
from ete3 import NCBITaxa
import six


class CorrInfo:
    def __init__(self, line, metabolite, cluster, pval, corr, pair):
        self.metabolite = metabolite
        self.cluster = cluster
        self.line = line
        self.pval = pval
        self.corr = corr
        self.pair = pair
        self.ctrl = None


def read_corr(corr_file):
    result = {}
    with open(corr_file) as f:
        for l in f:
            line = l.strip()
            if line:
                if "pval" not in line:
                    data = line.split(";")
                    metabolite = data[0]
                    cluster_no = data[1].split("_")[-1]
                    if cluster_no not in result:
                        result[cluster_no] = {}
                    if metabolite not in result[cluster_no]:
                        result[cluster_no][metabolite] = set()
                    result[cluster_no][metabolite].add(
                        CorrInfo(line=line.strip(), metabolite=data[0], cluster=data[1], pval=float(data[3]),
                                 corr=float(data[3]), pair=data[2]))
    return result


def compare_cntrl_corr(corr_file, cntrl_corr_info):
    for cluster_no, metabolit_relation in corr_file.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            for cntrl_data in cntrl_corr_info[cluster_no][metabolite]:
                for searched_data in metabolite_cl_data:
                    if cntrl_data.pair == searched_data.pair:
                        searched_data.ctrl = cntrl_data
    return corr_file


def get_database_sequences_info(fasta_database):
    result = {}
    with open(fasta_database) as f:
        for l in f:
            line = l.strip()
            if line.startswith(">"):
                uniprot_id = line.split("|")[1]
                protein_name = line.split("OS=")[0].split(" ", 1)[-1]
                organism_name = line.split("OS=", 1)[1].split("OX=")[0]
                organism_taxid = line.split("OX=")[-1].split(" ")[0]
                result[uniprot_id] = {"protein_name": protein_name, "organism_name": organism_name,
                                      "organism_taxid": organism_taxid}
    return result


def read_blast_result(blast_table, description, database_fasta_info, ncbi, result=None):
    # 5857992;4123967_2       tr|E3ZSC8|E3ZSC8_LISSE  57.143  112     48      0       1       112     1       112     2.67e-45        145
    if result is None:
        result = {}

    with open(blast_table) as f:
        for l in f:
            line = l.strip()
            if line:
                line = line.split()
                cl_no = line[0].split(";")[0]
                protein = line[1].split("|")[1]
                pident = float(line[2])
                print(database_fasta_info[protein]['organism_taxid'])
                lineage = ncbi.get_lineage(database_fasta_info[protein]['organism_taxid'])
                if 2759 not in lineage:
                    if cl_no in result:
                        if result[cl_no]["pident"] < pident:
                            result[cl_no]["pident"] = pident
                            result[cl_no]["protein"] = protein
                            result[cl_no]["description"] = description
                    else:
                        result[cl_no] = {"pident": pident, "protein": protein}
    return result


def fix_corr(corr_info, blast_result, database_info, save_old_line=True):
    for cluster_no, metabolit_relation in corr_info.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            for corr_data in metabolite_cl_data:
                if not save_old_line:
                    line = f"{corr_data.metabolite};{corr_data.cluster};{corr_data.pair};{corr_data.pval};{corr_data.corr}"
                    corr_data.line = line
                if cluster_no in blast_result:
                    corr_data.line += f";{blast_result[cluster_no]['pident']}"
                    corr_data.line += f";{blast_result[cluster_no]['protein']}"
                    corr_data.line += f";{database_info[blast_result[cluster_no]['protein']]['protein_name']}"
                    corr_data.line += f";{database_info[blast_result[cluster_no]['protein']]['organism_name']}"
                    if corr_data.ctrl is not None:
                        corr_data.line += f";{corr_data.ctrl.pval}"
                        corr_data.line += f";{corr_data.ctrl.corr}"
                else:
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    if corr_data.ctrl is not None:
                        corr_data.line += f";{corr_data.ctrl.pval}"
                        corr_data.line += f";{corr_data.ctrl.corr}"
    return corr_info


def save_corr(fixed_corr, out_file):
    with open(out_file, "w") as f:
        for cluset_no, metabolit_relation in fixed_corr.items():
            for metabolite, metabolite_cl_data in metabolit_relation.items():
                for corr_data in metabolite_cl_data:
                    f.write(corr_data.line)
                    f.write("\n")


@click.command()
@click.option('--corr_file', default="./", help='Folder with BLAST files.')
@click.option('--corr_ctrl_file', default="./", help='Folder with BLAST files.')
@click.option('--blast_files', default={}, help='')
@click.option('--fasta_database', default="./", help='')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
def main(corr_file, corr_ctrl_file, blast_files, fasta_database, out_file):
    blast_files = eval(blast_files)
    corr_info = read_corr(corr_file)
    import os
    if os.path.exists(corr_file):
        cntrl_corr_info = read_corr(corr_ctrl_file)
        corr_info = compare_cntrl_corr(corr_info, cntrl_corr_info)
    blast_result = {}
    database_fasta_info = get_database_sequences_info(fasta_database)
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    for blast_table, description in blast_files.items():
        if blast_table:
            blast_result = read_blast_result(blast_table, description, database_fasta_info, ncbi, blast_result)
    fixed_corr = fix_corr(corr_info, blast_result, database_fasta_info)
    save_corr(fixed_corr, out_file)


if __name__ == '__main__':
    main()
