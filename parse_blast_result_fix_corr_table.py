import os

import click
import numpy as np
from ete3 import NCBITaxa
import six


class CorrInfo:
    def __init__(self, line, metabolite, cluster, pval, corr, pair, metabo_values_list, gene_values_list):
        #             f.write(','.join([str(i) for i in datasets_data['metabo_values']]))
        #             f.write(";")
        #             f.write(','.join([str(i) for i in datasets_data['gene_values']]))
        self.metabolite = metabolite
        self.cluster = cluster
        self.line = line
        self.pval = pval
        self.corr = corr
        self.pair = pair
        self.metabo_values_list = metabo_values_list.split(',')
        self.gene_values_list = gene_values_list.split(',')
        self.ctrl = None
        self.fisher = None
        self.zou = None

    def count(self):
        fisher = independent_corr(self.corr, self.ctrl.corr, len(self.metabo_values_list),
                                  len(self.ctrl.metabo_values_list), method='fisher')
        zou = independent_corr(self.corr, self.ctrl.corr, len(self.metabo_values_list),
                               len(self.ctrl.metabo_values_list), method='zou')
        self.fisher = f"{fisher[0]};{fisher[1]}"
        self.zou = f"{zou[0]};{zou[1]}"


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
                        CorrInfo(line=line.strip(), metabolite=data[0], cluster=data[1], pval=float(data[2]),
                                 corr=float(data[3]), pair=data[2], metabo_values_list=data[4],
                                 gene_values_list=data[5]))
    return result


def compare_cntrl_corr(corr_file, cntrl_corr_info):
    for cluster_no, metabolit_relation in corr_file.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            for cntrl_data in cntrl_corr_info[cluster_no][metabolite]:
                for searched_data in metabolite_cl_data:
                    searched_data.ctrl = cntrl_data
                    searched_data.count()
    return corr_file


def get_database_sequences_info(fasta_database, result={}):
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


def read_blast_result(blast_table, description, database_fasta_info, ncbi, result=None, main_taxids=None):
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
                eval = float(line[10])
                # print(database_fasta_info[protein]['organism_taxid'])
                if main_taxids:
                    if database_fasta_info[protein]['organism_taxid'] in main_taxids:
                        if cl_no in result:
                            if result[cl_no]["pident"] < pident:
                                result[cl_no]["pident"] = pident
                                result[cl_no]["protein"] = protein
                                result[cl_no]["description"] = description
                                result[cl_no]["blast_table"] = blast_table
                                result[cl_no]["eval"] = eval
                            if result[cl_no]["pident"] == pident:
                                if eval < result[cl_no]["eval"]:
                                    result[cl_no]["pident"] = pident
                                    result[cl_no]["protein"] = protein
                                    result[cl_no]["description"] = description
                                    result[cl_no]["blast_table"] = blast_table
                                    result[cl_no]["eval"] = eval
                else:
                    try:
                        lineage = ncbi.get_lineage(database_fasta_info[protein]['organism_taxid'])
                        if 2759 not in lineage:
                            if cl_no in result:
                                if result[cl_no]["pident"] < pident:
                                    result[cl_no]["pident"] = pident
                                    result[cl_no]["protein"] = protein
                                    result[cl_no]["description"] = description
                                    result[cl_no]["blast_table"] = blast_table
                                    result[cl_no]["eval"] = eval
                                if result[cl_no]["pident"] == pident:
                                    if eval < result[cl_no]["eval"]:
                                        result[cl_no]["pident"] = pident
                                        result[cl_no]["protein"] = protein
                                        result[cl_no]["description"] = description
                                        result[cl_no]["blast_table"] = blast_table
                                        result[cl_no]["eval"] = eval
                            else:
                                result[cl_no] = {"pident": pident, "protein": protein, "description": description,
                                                 "blast_table": blast_table, "eval": eval}
                    except:
                        if cl_no in result:
                            if result[cl_no]["pident"] < pident:
                                result[cl_no]["pident"] = pident
                                result[cl_no]["protein"] = protein
                                result[cl_no]["description"] = description
                                result[cl_no]["blast_table"] = blast_table
                        else:
                            result[cl_no] = {"pident": pident, "protein": protein, "description": description,
                                             "blast_table": blast_table}
    return result


import numpy as np
from scipy.stats import t, norm
from math import atanh, pow
from numpy import tanh


def rz_ci(r, n, conf_level=0.95):
    zr_se = pow(1 / (n - 3), .5)
    moe = norm.ppf(1 - (1 - conf_level) / float(2)) * zr_se
    zu = atanh(r) + moe
    zl = atanh(r) - moe
    return tanh((zl, zu))


def independent_corr(xy, ab, n, n2=None, twotailed=True, conf_level=0.95, method='fisher'):
    """
    Calculates the statistic significance between two independent correlation coefficients
    @param xy: correlation coefficient between x and y
    @param xz: correlation coefficient between a and b
    @param n: number of elements in xy
    @param n2: number of elements in ab (if distinct from n)
    @param twotailed: whether to calculate a one or two tailed test, only works for 'fisher' method
    @param conf_level: confidence level, only works for 'zou' method
    @param method: defines the method uses, 'fisher' or 'zou'
    @return: z and p-val
    """

    if method == 'fisher':
        xy_z = 0.5 * np.log((1 + xy) / (1 - xy))
        ab_z = 0.5 * np.log((1 + ab) / (1 - ab))
        if n2 is None:
            n2 = n

        se_diff_r = np.sqrt(1 / (n - 3) + 1 / (n2 - 3))
        diff = xy_z - ab_z
        z = abs(diff / se_diff_r)
        p = (1 - norm.cdf(z))
        if twotailed:
            p *= 2

        return z, p
    elif method == 'zou':
        L1 = rz_ci(xy, n, conf_level=conf_level)[0]
        U1 = rz_ci(xy, n, conf_level=conf_level)[1]
        L2 = rz_ci(ab, n2, conf_level=conf_level)[0]
        U2 = rz_ci(ab, n2, conf_level=conf_level)[1]
        lower = xy - ab - pow((pow((xy - L1), 2) + pow((U2 - ab), 2)), 0.5)
        upper = xy - ab + pow((pow((U1 - xy), 2) + pow((ab - L2), 2)), 0.5)
        return lower, upper


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
                    corr_data.line += f";https://www.uniprot.org/uniprotkb/{blast_result[cluster_no]['protein']}/entry"
                    corr_data.line += f";{blast_result[cluster_no]['description']}"
                    corr_data.line += f";{blast_result[cluster_no]['blast_table']}"
                    corr_data.line += f";{database_info[blast_result[cluster_no]['protein']]['protein_name']}"
                    corr_data.line += f";{database_info[blast_result[cluster_no]['protein']]['organism_name']}"
                    corr_data.line += f";{database_info[blast_result[cluster_no]['protein']]['organism_taxid']}"

                    if corr_data.ctrl is not None:
                        corr_data.line += f";{corr_data.ctrl.pval}"
                        corr_data.line += f";{corr_data.ctrl.corr}"
                        corr_data.line += f";{corr_data.zou}"
                        corr_data.line += f";{corr_data.fisher}"
                else:
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    corr_data.line += f";"
                    if corr_data.ctrl is not None:
                        corr_data.line += f";{corr_data.ctrl.pval}"
                        corr_data.line += f";{corr_data.ctrl.corr}"
                        corr_data.line += f";{corr_data.zou}"
                        corr_data.line += f";{corr_data.fisher}"
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
@click.option('--fasta_databases', default="./", help='')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
@click.option('--main_taxids', default=None, help='Out file with genes statistics.')
def main(corr_file, corr_ctrl_file, blast_files, fasta_databases, out_file, main_taxids):
    blast_files = eval(blast_files)
    corr_info = read_corr(corr_file)
    import os
    if os.path.exists(corr_ctrl_file):
        cntrl_corr_info = read_corr(corr_ctrl_file)
        corr_info = compare_cntrl_corr(corr_info, cntrl_corr_info)
    blast_result = {}
    database_fasta_info = {}
    for fasta_db in fasta_databases.split(","):
        database_fasta_info = get_database_sequences_info(fasta_db, database_fasta_info)
    ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()
    for blast_table, description in blast_files.items():
        if blast_table:
            blast_result = read_blast_result(blast_table, description, database_fasta_info, ncbi, blast_result,
                                             main_taxids)
    fixed_corr = fix_corr(corr_info, blast_result, database_fasta_info)
    save_corr(fixed_corr, out_file)


if __name__ == '__main__':
    main()
