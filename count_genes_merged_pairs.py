import os

import click
import numpy as np


@click.command()
@click.option('--in_folder', default="./", help='Folder with BLAST files.')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
@click.option('--preffix_files', default="./", help='Preffix of files to analyse.')
def main(out_file, in_folder, preffix_files):
    file_list = [i for i in os.listdir(in_folder) if i.startswith(preffix_files)]
    res = {}
    genes = set()
    print("start")
    #print(file_list)
    with open(out_file, "w") as f: 
        f.write("genes;")
        f.write(";".join(file_list))
        f.write("\n")
        for e,file_name in enumerate(file_list):
            print(e, file_name)
            probe_id=file_name.split(".")
            res[probe_id] = {}
            path_to_csv = os.path.join(in_folder, file_name)
            data = np.loadtxt(path_to_csv, dtype=str)
            for line in data:
                if float(line[2]) > 90 and float(line[10]) < 0.000001:
                    if line[1] not in res[probe_id]:
                        res[probe_id][line[1]] = {line[0]}
                    else:
                        res[probe_id][line[1]].add(line[0])
                    genes.add(line[1])
            res[probe_id]={i:len(j) for i,j in res[file_name].items()}
        for gene in list(genes):
            f.write(f"{str(gene)};")
            for file_name in file_list:
                f.write(f"{str(res[probe_id].get(gene, 0))};")
            f.write("\n")

if __name__=='__main__':
    main()
