import os

import click
import numpy as np

def read_corr(corr_file):
    result={}
    with open(corr_file) as f:
        for l in f:
            line=l.strip()
            if line:
                data=line.split(";")
                metabolite=data[0]
                cluster_no=data[1].split("_")[-1]
                if cluster_no not in result:
                    result[cluster_no]={}
                if metabolite not in result[cluster_no]:
                    result[cluster_no][metabolite]={}
                result[cluster_no][metabolite]["line"]=line.strip()
                result[cluster_no][metabolite]["pval"]=float(data[4])
                result[cluster_no][metabolite]["corr"]=float(data[3])

    return result
                

def compare_cntrl_corr(corr_file, cntrl_corr_info):
    for cluset_no, metabolit_relation in corr_file.items():
        for metabolite, metabolite_cl_data in metabolit_relation.item():
            metabolite_cl_data["pval_cntrl"]=cntrl_corr_info[cluster_no][metabolite]["pval"]
            metabolite_cl_data["corr_cntrl"]=cntrl_corr_info[cluster_no][metabolite]["corr"]
    return corr_file

def get_database_sequences_info(fasta_database):
    result={}
    with open(fasta_database) as f:
        for l in f:
            line=l.strip()
            if line.startswith(">"):
                uniprot_id=line.split("|")[1]
                protein_name=line.split("OS=")[0].split(" ", 1)[-1]
                organism_name=line.split("OS=", 1)[1].split("OX=")[0]
                result[uniprot_id]={"protein_name":protein_name, "organism_name":organism_name}
    return result
                
def read_blast_result(blast_table):
    # 5857992;4123967_2       tr|E3ZSC8|E3ZSC8_LISSE  57.143  112     48      0       1       112     1       112     2.67e-45        145
    result={}
    with open(blast_table) as f:
        for l in f:
            line=l.strip()
            if line:
                line=line.split()
                cl_no=line[0].split(";")[0]
                protein=line[1].split("|")[1]
                pident=float(line[2])
                if cl_no in result:
                    if result[cl_no]["pident"]>pident:
                        result[cl_no]["pident"]=pident
                        result[cl_no]["protein"]=protein
                else:
                    result[cl_no]={"pident":pident, "protein":protein}
    return result
    
def fix_corr(corr_info, database_info):
    for cluster_no, metabolit_relation in corr_info.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            metabolite_cl_data["line"]+=f";{database_info[cluster_no]['pident']}"
            metabolite_cl_data["line"]+=f";{database_info[cluster_no]['protein']}"
            metabolite_cl_data["line"]+=f";{metabolite_cl_data['pval_cntrl']}"
            metabolite_cl_data["line"]+=f";{metabolite_cl_data['corr_cntrl']}"
    return corr_file
    
def save_corr(fixed_corr, out_file):
    with open(out_file, "w") as f:
        for cluset_no, metabolit_relation in fixed_corr.items():
            for metabolite, metabolite_cl_data in metabolit_relation.item():
                f.write(metabolite_cl_data["line"])
                f.write("\n")


@click.command()
@click.option('--corr_file', default="./", help='Folder with BLAST files.')
@click.option('--corr_ctrl_file', default="./", help='Folder with BLAST files.')
@click.option('--blast_table', default="./", help='')
@click.option('--fasta_database', default="./", help='')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
def main(corr_file, corr_ctrl_file, blast_table, fasta_database, out_file):
    corr_info = read_corr(corr_file)
    cntrl_corr_info = read_corr(corr_ctrl_file)
    blast_result = read_blast_result(blast_table)
    database_fasta_info = get_database_sequences_info(fasta_database)
    fixed_corr = fix_corr(corr_info, database_fasta_info)
    save_corr(fixed_corr, out_file)


if __name__=='__main__':
    main()
