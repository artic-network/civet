#!/usr/bin/env python3
from civet.utils import log_colours as colour
from civet.utils import misc
import sys
import os
import csv
from Bio import SeqIO


def check_dir_for_file(file_list, extension,argument,key,config):
    match = []
    for fn in file_list:
        if fn.endswith(extension):
            match.append(fn)
    if len(match) > 1:
        sys.stderr.write(colour.cyan(f"Error: More than 1 {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
        sys.exit(-1)
    elif len(match) == 0:
        sys.stderr.write(colour.cyan(f"Error: No {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
        sys.exit(-1)
    elif len(match) == 1:
        config[key] = match[0]


def background_data_load(config):
    
    file_list = []
    for r,d,f in os.walk(config["datadir"]):
        for fn in f:
            file_list.append(os.path.join(r,fn))

    if not "background_csv" in config:
        check_dir_for_file(file_list, "csv","-bc/--background-csv","background_csv",config)
    if not "background_fasta" in config:
        check_dir_for_file(file_list, "fasta","-bf/--background-fasta","background_fasta",config)


def data_group_parsing(datadir,background_csv,background_fasta,data_column,config):
    
    if datadir:
        datadir = os.path.abspath(datadir)
    #try climb datadir

    misc.add_arg_to_config("datadir",datadir,config)
    misc.add_arg_to_config("data_column",data_column,config)

    misc.add_file_to_config("background_csv",background_csv,config)
    misc.add_file_to_config("background_fasta",background_fasta,config)

    background_data_load(config)

