#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import sys
import os
import csv
from Bio import SeqIO

def check_coords_within_reference_length(config):
    
    reference_fasta = config["reference_fasta"]
    ref_length = 0
    num_records = 0

    ts = config["trim_start"]
    te = config["trim_end"]

    try:
        ts = int(ts)
        config["trim_start"] = ts
    except:
        sys.stderr.write(cyan(f"`-ts/--trim-start` must be an integer.\n"))
        sys.exit(-1)
    
    try:
        te = int(te)
        config["trim_end"] = te
    except:
        sys.stderr.write(cyan(f"`-te/--trim-end` must be an integer.\n"))
        sys.exit(-1)

    try:
        for record in SeqIO.parse(reference_fasta,"fasta"):
            num_records +=1
            ref_length = len(record)
        if num_records > 1:
            sys.stderr.write(cyan(f"Please provide one sequence in reference fasta file: ") + f"{reference_fasta}\n")
            sys.exit(-1)
        if num_records == 0:
            sys.stderr.write(cyan(f"No sequences found in reference fasta file: ") + f"{reference_fasta}\n")
            sys.exit(-1)
    except:
        sys.stderr.write(cyan(f"Cannot read reference fasta file:\n") + f"{reference_fasta}\n" + cyan("Please check file is in correct format.\n"))
        sys.exit(-1)
    
    if int(te) < int(ts):
        sys.stderr.write(cyan(f"`-te/--trim-end` coordinate must be greater than `-ts/--trim-start` coordinate.\n"))
        sys.exit(-1)

    if int(te) <= ref_length and int(ts) <= ref_length:
        pass
    else:
        sys.stderr.write(cyan(f"`-te/--trim-end` and `-ts/--trim-start` must be within the length of the reference genome:\n") + f"Genome: {reference_fasta}\nTrim-start: {te}\nTrim-end: {ts}\n" + cyan("Please check file is in correct format.\n"))
        sys.exit(-1)

def check_catchment_configuration(config):

    try:
        cs = int(config["catchment_size"])
        config["catchment_size"] = cs
    except:
        sys.stderr.write(cyan(f"`-cs/--catchment_size` must be an integer.\n"))
        sys.exit(-1)

def analysis_group_parsing(reference_fasta,trim_start,trim_end,catchment_size,downsample,query_limit,config):
    """
    parses the data group arguments 
    --datadir (Default $DATADIR)
    --background-csv (Default False)
    --background-fasta (Default False)
    --data-column (Default central_sample_id)
    """

    # if command line arg, overwrite config value
    misc.add_arg_to_config("trim_start",trim_start,config)
    misc.add_arg_to_config("trim_end",trim_end,config)
    misc.add_file_to_config("reference_fasta",reference_fasta,config)
    misc.add_arg_to_config("catchment_size",catchment_size,config)
    misc.add_arg_to_config("downsample",downsample,config)
    misc.add_arg_to_config("query_limit",query_limit,config)

    check_coords_within_reference_length(config)
    check_catchment_configuration(config)
