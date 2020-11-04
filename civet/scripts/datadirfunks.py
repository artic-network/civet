#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import date
from collections import defaultdict
import pandas as pd
import random

import tempfile
import pkg_resources
import yaml

from reportfunk.funks import io_functions as qcfunk
from reportfunk.funks import prep_data_functions as prep_data
from reportfunk.funks import table_functions as table_func

today = date.today()

def get_background_files(data_dir,background_metadata):
    background_seqs = ""
    background_tree = ""
    data_date = ""
    for r,d,f in os.walk(data_dir):
        for fn in f:
            if fn.endswith(".fasta"):
                background_seqs = os.path.join(data_dir, fn)
                data_date = fn.split("_")[2]
                if not data_date.startswith("20"):
                    data_date = ""
            elif fn.endswith(".nexus"):
                background_tree = os.path.join(data_dir, fn)
            elif background_metadata == "" and fn.endswith(".csv"):
                background_metadata = os.path.join(data_dir, fn)

    return background_seqs, background_tree, background_metadata, data_date
    

def print_data_error(data_dir):
    sys.stderr.write(qcfunk.cyan(f"Error: data directory should contain the following files or additionally supply a background metadata file:\n") + f"\
    - global_2020-XX-YY_tree.nexus\n\
    - global_2020-XX-YY_metadata.csv\n\
    - global_2020-XX-YY_alignment.fasta\n"+qcfunk.cyan(f"\
To run civet please specify a local directory with the appropriate files, optionally supply a custom metadata file\n\n"""))


def get_datadir(args_datadir,args_metadata,cwd,config):
    data_dir = ""
    background_metadata = ""

    if args_metadata:
        expanded_path = os.path.expanduser(args_metadata)
        background_metadata = os.path.join(cwd, expanded_path)
        if not os.path.exists(background_metadata):
            sys.stderr.write(qcfunk.cyan(f"Error: can't find metadata file at {background_metadata}.\n"))
            sys.exit(-1)

    elif "background_metadata" in config:
        if config["background_metadata"]:
            expanded_path = os.path.expanduser(config["background_metadata"])
            background_metadata = os.path.join(config["path_to_query"], expanded_path)
            if not os.path.exists(background_metadata):
                sys.stderr.write(qcfunk.cyan(f"Error: can't find metadata file at {background_metadata}.\n"))
                sys.exit(-1)
            
    
    elif args_datadir:
        data_dir = os.path.join(cwd, args_datadir)

    elif "datadir" in config:
        if config["datadir"]:
            expanded_path = os.path.expanduser(config["datadir"])
            data_dir = os.path.join(config["path_to_query"], expanded_path)
        else:
            data_dir = os.path.join(cwd, "civet-cat")

    if not os.path.exists(data_dir):
        print_data_error(data_dir)
        sys.exit(-1)
        
    background_seqs, background_tree, background_metadata, data_date = get_background_files(data_dir,background_metadata)

    config["datadir"] = data_dir
    config["data_date"] = data_date

    if not os.path.isfile(background_tree) or not os.path.isfile(background_seqs) or not os.path.isfile(background_metadata):
        print_data_error(data_dir)
        sys.exit(-1)
    else:
        config["background_metadata"] = background_metadata
        config["background_seqs"] = background_seqs
        config["background_tree"] = background_tree

        print("Found data:")
        print("    -",background_seqs)
        print("    -",background_metadata)
        print("    -",background_tree,"\n")
