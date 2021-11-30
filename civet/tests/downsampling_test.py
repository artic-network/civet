import collections
import csv
from Bio import SeqIO
import hashlib
from numpy.random import choice
import sys
from civet.utils.log_colours import green,cyan,red
from datetime import datetime
import os
import yaml

from civet import catchment_parsing

def test_downsampling(configfile):
    print(datetime.now())
    with open(configfile, 'r') as f:
        config_loaded = yaml.safe_load(f)

    fasta  = os.path.join(config_loaded["tempdir"],"hashed.aln.fasta")
    in_csv = os.path.join(config_loaded["data_outdir"],"catchments","query_metadata.catchments.csv")
    catchment_dir = os.path.join(config_loaded["data_outdir"],"catchments")
    out_csv = os.path.join(config_loaded["data_outdir"],"catchments","master_metadata.downsample.csv")
    print(datetime.now())
    catchment_parsing.downsample_if_building_trees(in_csv, out_csv, config_loaded)
    print(datetime.now())
    if '3' in config_loaded["report_content"]:
        print(green("Writing catchment fasta files."))
        catchment_parsing.write_catchment_fasta(out_csv,fasta,catchment_dir,config_loaded)
        print(datetime.now())

test_downsampling("/Users/s1680070/repositories/civet3/civet_2021-06-29/catchments/config.yaml")