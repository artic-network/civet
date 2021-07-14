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
        mq = int(config["max_queries"])
        config["max_queries"] = mq
    except:
        sys.stderr.write(cyan(f"`-mq/--max-queries` must be an integer.\n"))
        sys.exit(-1)

def check_for_background_header(config):
    with open(config["background_metadata"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if not config["downsample_column"] in header:
            sys.stderr.write(cyan(f"`--ds/--downsample` column specified ({config['downsample_column']}), but is not found in background metadata file. Please indicate a valid metadata column.\n"))
            sys.exit(-1)

def check_max_memory(config):
    try:
        mm = int(config["max_memory"])
        config["max_memory"] = mm
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
        config["mode"]="random"
        print(green("Downsampling strategy: ") + "random")

    elif mode=="normalise":
        config["mode"]="normalise"

        if not col:
            sys.stderr.write(cyan(f"`--ds/--downsample` mode `normalise` specified, but no column indicated. Please indicate which metadata column with `--downsample mode=normalise <column>`, where <column> is the column you want to normalise by.\n"))
            sys.exit(-1)
        else:
            config["downsample_column"]=col

        print(green("Downsampling strategy: ") + "normalise")
        print(green("Normalising over: ") + col)

    elif mode=="enrich":
        config["mode"]="enrich"
        if factor:
            try:
                factor = float(factor)
            except:
                sys.stderr.write(cyan(f"`--ds/--downsample` factor must be numerical.\n"))
                sys.exit(-1)
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

        print(green("Downsampling strategy: ") + "enrich")
        print(green("Enriching for ") + downsample_column + green(" with ") + field + green(" by a factor of ") + str(config['factor']) + ".")
    else:
        sys.stderr.write(cyan(f"`--ds/--downsample` mode not one of random, normalise or enrich.\n"))
        sys.exit(-1)
    
    if config["downsample_column"]:
        check_for_background_header(config)


def analysis_group_parsing(reference_fasta,trim_start,trim_end,max_queries,max_memory,config):

    # if command line arg, overwrite config value
    misc.add_arg_to_config("trim_start",trim_start,config)
    misc.add_arg_to_config("trim_end",trim_end,config)
    misc.add_file_to_config("reference_fasta",reference_fasta,config)
    misc.add_arg_to_config("max_queries",max_queries,config)
    misc.add_arg_to_config("max_memory",max_memory,config)

    check_coords_within_reference_length(config)
    check_query_limit(config)
    check_max_memory(config)

def check_if_int(key,config):
    if config[key]:
        try:
            config[key] = int(config[key])
        except:
            sys.stderr.write(cyan(f"`{key}` must be numerical.\n"))
            sys.exit(-1)

def distance_configuration(config):
    if config["snp_distance"]:
        check_if_int("snp_distance",config)
        print(green("Overwriting SNP radius for catchment with distance: ") + f"{config['snp_distance']}")
        config["snp_distance_up"] = config["snp_distance"]
        config["snp_distance_down"] = config["snp_distance"]
        config["snp_distance_side"] = config["snp_distance"]
    else:
        check_if_int("snp_distance_up",config)
        check_if_int("snp_distance_down",config)
        check_if_int("snp_distance_side",config)
        print(green("SNP distance radius for catchment:") + f"\n\t- Up {config['snp_distance_up']}\n\t- Down {config['snp_distance_down']}\n\t- Side {config['snp_distance_side']}")

def catchment_group_parsing(catchment_size,downsample,distance,distance_up,distance_down,distance_side,config):

    misc.add_arg_to_config("catchment_size",catchment_size,config)
    misc.add_arg_to_config("downsample",downsample,config)

    misc.add_arg_to_config("snp_distance",distance,config)
    misc.add_arg_to_config("snp_distance_up",distance_up,config)
    misc.add_arg_to_config("snp_distance_down",distance_down,config)
    misc.add_arg_to_config("snp_distance_side",distance_side,config)

    check_catchment_configuration(config)
    parse_downsampling_config(config)
    distance_configuration(config)