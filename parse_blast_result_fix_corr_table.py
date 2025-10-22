import os

import click
import numpy as np


@click.command()
@click.option('--in_folder', default="./", help='Folder with BLAST files.')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
@click.option('--preffix_files', default="./", help='Preffix of files to analyse.')
def main():
    pass


if __name__=='__main__':
    main()
