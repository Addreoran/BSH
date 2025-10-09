import os

import click
import numpy as np


@click.command()
@click.option('--in_folder', default=1, help='Number of greetings.')
@click.option('--out_file', default=1, help='Number of greetings.')
@click.option('--preffix_files', default=1, help='Number of greetings.')
def main(out_file, in_folder, preffix_files):
    file_list = os.listdir(in_folder)
    res = {}
    genes = set()

    for file_name in file_list:
        if file_name.startswith(preffix_files):
            res[file_name] = {}
            path_to_csv = os.path.join(in_folder, file_name)
            data = np.loadtext(path_to_csv, dtype=str, names=False)
            for line in data:
                if float(line[2]) > 90 and float(line[10]) < 0.000001:
                    if line[1] not in res[file_name]:
                        res[file_name][line[1]] = {line[0]}
                    else:
                        res[file_name][line[1]].add(line[0])
                    genes.add(line[0])

    with open(out_file, "w") as f:
        f.write("genes;")
        f.write(";".join(file_list))
        f.write("\n")
        for gene in list(genes):
            f.write(f"{str(gene)};")
            for file_name in file_list:
                f.write(f"{str(res[file_name].get(gene, 0))};")
            f.write("\n")
