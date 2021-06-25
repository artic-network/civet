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
        sys.stderr.write(cyan(f"`-cs/--catchment-size` must be an integer.\n"))
        sys.exit(-1)

def check_query_limit(config):
    try:
        ql = int(config["query_limit"])
        config["query_limit"] = cs
    except:
        sys.stderr.write(cyan(f"`-ql/--query-limit` must be an integer.\n"))
        sys.exit(-1)

def check_for_background_header(config):
    with open(config["background_csv"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if not config["downsample_column"] in header:
            sys.stderr.write(cyan(f"`--ds/--downsample` column specified ({config['downsample_column']}), but is not found in background metadata file. Please indicate a valid metadata column.\n"))
            sys.exit(-1)

def check_max_memory(config):
    try:
        ql = int(config["max_memory"])
        config["max_memory"] = cs
    except:
        sys.stderr.write(cyan(f"`-mem/--max-memory` must be an integer.\n"))
        sys.exit(-1)

def parse_downsampling_config(config):
    downsample = config["downsample"]

    if not type(downsample)== list:
        downsample = downsample.split(" ")
    
    if len(downsample)>3:
        sys.stderr.write(cyan(f"`--ds/--downsample` supplied with more arguments than expected.\n"))
        sys.exit(-1)

    mode = ''
    factor = ''
    col = ''
    for i in downsample:
        if i.startswith("mode"):
            mode = i.split("=")[1]
        elif i.startswith("factor"):
            factor = i.split("=")[1]
        else:
            col = i

    if not mode:
        sys.stderr.write(cyan(f"`mode` not indicated for `--ds/--downsample`. Please indicate one of random, normalise or enrich with `--downsample mode=<random,normalise,enrich>`\n"))
        sys.exit(-1)
    
    elif mode=="random":
        config["mode"]=="random"

    elif mode=="normalise":
        config["mode"]=="normalise"

        if not col:
            sys.stderr.write(cyan(f"`--ds/--downsample` mode `normalise` specified, but no column indicated. Please indicate which metadata column with `--downsample mode=normalise <column>`, where <column> is the column you want to normalise by.\n"))
            sys.exit(-1)
        else:
            config["downsample_column"]=col

    elif mode=="enrich":
        config["mode"]=="enrich"
        if factor:
            config["factor"]=factor
        else:
            config["factor"]=10
        if not col:
            sys.stderr.write(cyan(f"`--ds/--downsample` mode `enrich` specified, but no column or field indicated. Please indicate which metadata column and field with <column>=<field>, for example `--downsample mode=enrich factor=10 country=UK`, where country is the column to normalise by and UK is the field you'd like enriched for.\n"))
            sys.exit(-1)
        elif '=' not in col:
            sys.stderr.write(cyan(f"`--ds/--downsample` mode `enrich` specified, but no field indicated. Please indicate which metadata column and field with <column>=<field>, for example: `--downsample mode=enrich factor=10 country=UK`, where country is the column to normalise by and UK is the field you'd like enriched for.\n"))
            sys.exit(-1)
        else:
            downsample_column,field=col.split("=")
            config["downsample_column"]=downsample_column
            config["downsample_field"]=field
    else:
        sys.stderr.write(cyan(f"`--ds/--downsample` mode not one of random, normalise or enrich.\n"))
        sys.exit(-1)
    
    if config["downsample_column"]:
        check_for_background_header(config)


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
    misc.add_arg_to_config("max_memory",max_memory,config)
    misc.add_arg_to_config("downsample",downsample,config)

    check_coords_within_reference_length(config)
    check_catchment_configuration(config)
    check_query_limit(config)
    check_max_memory(config)

    parse_downsampling_config(config)
