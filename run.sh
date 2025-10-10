python3 ./count_genes.py --in_folder /data02/maria/kwasy/shotgun_BSH/ --out_file ./result.csv --preffix_files mgshot_
python3 ./merge_by_cd_hit_res.py --out_file ./result_merged.csv --count_genes_file ./result.csv --cd_hit_res /data02/maria/kwasy/shotgun_BSH/clusters.clst
