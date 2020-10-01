#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import date
from collections import defaultdict
import pandas as pd

import tempfile
import pkg_resources
import yaml

from reportfunk.funks import io_functions as qcfunk
from reportfunk.funks import prep_data_functions as prep_data
from reportfunk.funks import table_functions as table_func

today = date.today()

def get_defaults():
    default_dict = {
                    "title": "# Cluster investigation",
                    "outbreak_id": "",
                    "report_date": today,# date investigation was opened
                    "authors": "", # List of authors, affiliations and contact details
                    "description": "",
                    "conclusions": "",
                    "max_ambiguity":0.5,
                    "min_length":10000,
                    "datadir":"civet-cat",
                    "input_column":"name",
                    "data_column":"central_sample_id",
                    "sample_date_column":"sample_date",
                    "display_name":False,
                    "private":True,
                    "distance":2,
                    "up_distance":False,
                    "down_distance":False,
                    "collapse_threshold":1,
                    "sequencing_centre":"DEFAULT",
                    "tree_fields":"adm1",
                    "local_lineages":False,
                    "map_sequences":False,
                    "map_info":False,
                    "input_crs":False,
                    "colour_map_by":False,
                    "date_restriction":False,
                    "date_range_start":False,
                    "date_range_end":False,
                    "node_summary":"country",
                    "date_window":7,
                    "colour_by":"adm1=viridis",
                    "label_fields":False,
                    "date_fields":False,
                    "no_snipit":False,
                    "include_snp_table":False,
                    "include_bars":False,
                    "cog_report":False,
                    "omit_appendix":False,
                    "table_fields":["sample_date", "uk_lineage", "lineage", "phylotype"],
                    "threads":1,
                    "force":True,
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674
                    }
    return default_dict

def define_seq_db(config,default_dict):
    config["seq_db"] = config["background_seqs"]
    

def get_package_data(thisdir,config):
    reference_fasta = pkg_resources.resource_filename('civet', 'data/reference.fasta')
    outgroup_fasta = pkg_resources.resource_filename('civet', 'data/outgroup.fasta')
    polytomy_figure = pkg_resources.resource_filename('civet', 'data/polytomies.png')
    report_args = pkg_resources.resource_filename('civet', 'data/report_arguments.txt')
    footer_fig = pkg_resources.resource_filename('civet', 'data/footer.png')
    clean_locs = pkg_resources.resource_filename('civet', 'data/mapping_files/adm2_cleaning.csv')
    map_input_1 = pkg_resources.resource_filename('civet', 'data/mapping_files/gadm36_GBR_2.json')
    map_input_2 = pkg_resources.resource_filename('civet', 'data/mapping_files/channel_islands.json')  
    map_input_3 = pkg_resources.resource_filename('civet', 'data/mapping_files/NI_counties.geojson')  
    map_input_4 = pkg_resources.resource_filename('civet', 'data/mapping_files/Mainland_HBs_gapclosed_mapshaped_d3.json')
    map_input_5 = pkg_resources.resource_filename('civet', 'data/mapping_files/urban_areas_UK.geojson')
    map_input_6 = pkg_resources.resource_filename('civet', 'data/mapping_files/UK_outPC_coords.csv')
    spatial_translations_1 = pkg_resources.resource_filename('civet', 'data/mapping_files/HB_Translation.pkl')
    spatial_translations_2 = pkg_resources.resource_filename('civet', 'data/mapping_files/adm2_regions_to_coords.csv')
    appendix_text = pkg_resources.resource_filename('civet', 'data/appendix.txt')
    config["reference_fasta"] = reference_fasta
    config["outgroup_fasta"] = outgroup_fasta
    config["polytomy_figure"] = polytomy_figure
    config["report_args"] = report_args
    config["footer"] = footer_fig
    config["appendix"] = appendix_text
    
    config["clean_locs"] = clean_locs
    config["uk_map"] = map_input_1
    config["channels_map"] = map_input_2
    config["ni_map"] = map_input_3
    config["uk_map_d3"] = map_input_4
    config["urban_centres"] = map_input_5
    config["pc_file"] = map_input_6
    config["HB_translations"] = spatial_translations_1
    config["PC_translations"] = spatial_translations_2

    report_template = os.path.join(thisdir, 'scripts','civet_template.pmd')

    if not os.path.exists(report_template):
        sys.stderr.write(qcfunk.cyan(f'Error: cannot find report_template at {report_template}\n'))
        sys.exit(-1)

    config["report_template"] = report_template

def print_data_error(data_dir):
    sys.stderr.write(qcfunk.cyan(f"Error: data directory should contain the following files or additionally supply a background metadata file:\n") + f"\
    - cog_global_tree.nexus\n\
    - cog_global_metadata.csv\n\
    - cog_global_alignment.fasta\n"+qcfunk.cyan(f"\
To run civet please either\n1) ssh into CLIMB and run with --CLIMB flag\n\
2) Run using `--remote-sync` flag and your CLIMB username specified e.g. `-uun climb-covid19-otoolexyz`\n\
3) Specify a local directory with the appropriate files, optionally supply a custom metadata file\n\n"""))

def rsync_data_from_climb(uun, data_dir):
    rsync_command = f"rsync -avzh {uun}@bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{data_dir}'"
    print(qcfunk.green(f"Syncing civet data to {data_dir}"))
    status = os.system(rsync_command)
    if status != 0:
        sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your user name is a valid CLIMB username e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and are in the UK\n\n"))
        sys.exit(-1)

def get_remote_data(uun,background_metadata,data_dir,config):
    config["remote"]= True
    head,tail = os.path.split(data_dir)
    if tail == "civet-cat":
        path_for_syncing = head
    else:
        path_for_syncing = data_dir
    if uun:
        config["username"] = uun
        rsync_data_from_climb(uun, path_for_syncing)
    elif "username" in config:
        uun = config["username"]
        rsync_data_from_climb(uun, path_for_syncing)
    elif "uun" in config:
        uun = config["uun"]
        rsync_data_from_climb(uun, path_for_syncing)
    else:
        rsync_command = f"rsync -avzh bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{path_for_syncing}'"
        print(f"Syncing civet data to {path_for_syncing}")
        status = os.system(rsync_command)
        if status != 0:
            sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your ssh is configured with Host bham.covid19.climb.ac.uk\nAlternatively enter your CLIMB username with -uun e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and check if you are in the UK\n\n"))
            sys.exit(-1)

    background_seqs = ""
    background_tree = ""
    
    if not background_metadata:
        background_metadata = os.path.join(path_for_syncing,"civet-cat","cog_global_metadata.csv")
    
    background_seqs= os.path.join(path_for_syncing,"civet-cat","cog_global_alignment.fasta")
    background_tree = os.path.join(path_for_syncing,"civet-cat","cog_global_tree.nexus")

    config["datadir"] = os.path.join(path_for_syncing,"civet-cat")
    if not os.path.exists(config["datadir"]):
        print(qcfunk.cyan(f"Error: data directory not found at {data_dir}.\n"))
        sys.exit(-1)

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

def get_datadir(args_climb,args_uun,args_datadir,args_metadata,remote,cwd,config,default_dict):
    data_dir = ""
    background_metadata = ""

    if args_metadata:
        background_metadata = os.path.join(cwd, args_metadata)
        if not os.path.exists(background_metadata):
            sys.stderr.write(qcfunk.cyan(f"Error: can't find metadata file at {background_metadata}.\n"))
            sys.exit(-1)

    elif "background_metadata" in config:
        expanded_path = os.path.expanduser(config["background_metadata"])
        background_metadata = os.path.join(config["path_to_query"], expanded_path)
        if not os.path.exists(background_metadata):
            sys.stderr.write(qcfunk.cyan(f"Error: can't find metadata file at {background_metadata}.\n"))
            sys.exit(-1)
            
    if args_climb:
        data_dir = "/cephfs/covid/bham/civet-cat"
        if os.path.exists(data_dir):
            config["remote"] = False
            config["username"] = ""
        else:
            sys.stderr.write(qcfunk.cyan(f"Error: --CLIMB argument called, but CLIMB data path doesn't exist.\n"))
            sys.exit(-1)

    elif args_datadir:
        data_dir = os.path.join(cwd, args_datadir)

    elif "datadir" in config:
        expanded_path = os.path.expanduser(config["datadir"])
        data_dir = os.path.join(config["path_to_query"], expanded_path)
    else:
        if remote:
            data_dir = cwd
        else:
            data_dir = os.path.join(cwd, default_dict["datadir"])

    if not remote:
        if not os.path.exists(data_dir):
            print_data_error(data_dir)
            sys.exit(-1)
            
        
        background_seqs = ""
        background_tree = ""
        
        if not background_metadata:
            background_metadata = os.path.join(data_dir,"cog_global_metadata.csv")
        background_seqs= os.path.join(data_dir,"cog_global_alignment.fasta")
        background_tree = os.path.join(data_dir,"cog_global_tree.nexus")

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

    elif remote:
        
        get_remote_data(args_uun, background_metadata, data_dir, config)

    config["datadir"]=data_dir


def prepping_civet_arguments(name_stem_input, tree_fields_input, graphic_dict_input, label_fields_input, date_fields_input, table_fields_input):

    tree_fields = prep_data.prep_argument_list(tree_fields_input)
    label_fields = prep_data.prep_argument_list(label_fields_input)
    date_fields = prep_data.prep_argument_list(date_fields_input)
    table_fields = prep_data.prep_argument_list(table_fields_input) 

    if "/" in name_stem_input:
        name_stem = name_stem_input.split("/")[-1]
    else:
        name_stem = name_stem_input

    graphic_dict = {}
    if type(graphic_dict_input) == str:
        splits = graphic_dict_input.split(",")
    else:
        splits = graphic_dict_input
    
    for element in splits:
        key = element.split(":")[0].replace(" ","").replace("'","")
        value = element.split(":")[1].replace(" ","").replace("'","")
        graphic_dict[key] = value
            
    for key in graphic_dict.keys():
        if key not in tree_fields:
            tree_fields.append(key)
  
    return name_stem, tree_fields, graphic_dict, label_fields, date_fields, table_fields

def local_lineages_qc(config,default_dict):

    query_file = config["query"]

    if config["local_lineages"]:

        if config["date_restriction"]:

            if config["date_range_start"] and type(config["date_range_start"]) == str:
                check_date_format(config["date_range_start"])
            if config["date_range_end"] and type(config["date_range_end"]) == str:
                check_date_format(config["date_range_end"])

            if config["date_range_start"] and config["date_range_end"]:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_range_start']} to {config['date_range_end']}"))
            elif config["date_range_start"]:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_range_start']} to present"))
            else:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_window']} days around the sampling range"))

def map_sequences_config(map_sequences,colour_map_by,map_inputs,input_crs,config,default_dict):
    
    if config["map_sequences"]:

        map_inputs = ""
        if config["map_info"]:
            map_inputs = config["map_info"].replace(" ","")
        else:
            sys.stderr.write(qcfunk.cyan('Error: coordinates or outer postcode not supplied for mapping sequences. Please provide either x and y columns as a comma separated string, or column header containing outer postcode.'))
            sys.exit(-1)

        if len(map_inputs.split(",")) == 2: #If x and y coordinates are provided
            if not config["input_crs"]:
                sys.stderr.write(qcfunk.cyan('Error: input coordinate system not provided for mapping. Please provide --input-crs eg EPSG:3395'))
                sys.exit(-1)
        else: #If an outer postcode column is provided
            config["input_crs"] = "EPSG:4326"
        
        with open(config["query"], newline="") as f:
            reader = csv.DictReader(f)
            column_names = reader.fieldnames

            relevant_cols = map_inputs.split(",")

            if config["colour_map_by"]:
                relevant_cols.append(colour_map_by)
            
            for map_arg in relevant_cols:
                map_arg = map_arg.replace(" ","")
                if map_arg not in column_names and map_arg not in config["background_metadata_headers"]:
                    sys.stderr.write(qcfunk.cyan(f"Error: {map_arg} column not found in metadata file or background database for mapping sequences"))
                    sys.exit(-1)

        if config["colour_map_by"]:
            if map_inputs == "adm2":
                print(qcfunk.cyan(f"NOTE: --colour-map-by not set up to colour by adm2. Please provide outer postcode or coordinates"))
            else:
                print(qcfunk.green(f"Colouring map by: " + f"{colour_map_by}"))
            config["colour_map_by"] = colour_map_by
            
        else:
            config["colour_map_by"] = False

        config["map_info"] = False
        config["input_crs"] = False
        config["colour_map_by"] = False


def get_sequencing_centre_header(config):
    
    sc_list = ["PHEC", 'LIVE', 'BIRM', 'PHWC', 'CAMB', 'NORW', 'GLAS', 'EDIN', 'SHEF',
                'EXET', 'NOTT', 'PORT', 'OXON', 'NORT', 'NIRE', 'GSTT', 'LOND', 'SANG',"NIRE"]

    sequencing_centre = config["sequencing_centre"]
    if sequencing_centre in sc_list or sequencing_centre == "DEFAULT":
        package_png = os.path.join("data","headers",f"{sequencing_centre}.png")
        sequencing_centre_source = pkg_resources.resource_filename('civet', package_png)
        print(qcfunk.green(f"Using header file from:") + f" {package_png}\n")
        config["sequencing_centre_source"] = sequencing_centre_source
        config["sequencing_centre_dest"] = os.path.join(config["outdir"],"figures",f"{sequencing_centre}.png")
        config["sequencing_centre_file"] = os.path.join(".","figures",f"{sequencing_centre}.png")
        config["sequencing_centre"] = sequencing_centre
    else:
        sc_string = "\n".join(sc_list)
        sys.stderr.write(qcfunk.cyan(f'Error: sequencing centre must be one of the following:\n{sc_string}\n'))
        sys.exit(-1)

def map_group_to_config(args,config,default_dict):

    ## local_lineages
    local_lineages = qcfunk.check_arg_config_default("local_lineages",args.local_lineages, config, default_dict)
    config["local_lineages"] = local_lineages

    ## date_restriction
    date_restriction = qcfunk.check_arg_config_default("date_restriction",args.date_restriction, config, default_dict)
    config["date_restriction"] = date_restriction

    ## date_range_start
    date_range_start = qcfunk.check_arg_config_default("date_range_start",args.date_range_start, config, default_dict)
    config["date_range_start"] = date_range_start

    ## date_range_end
    date_range_end = qcfunk.check_arg_config_default("date_range_end",args.date_range_end, config, default_dict)
    config["date_range_end"] = date_range_end

    ## date_window
    date_window = qcfunk.check_arg_config_default("date_window",args.date_window, config, default_dict)
    config["date_window"] = date_window

    ## map_sequences
    map_sequences = qcfunk.check_arg_config_default("map_sequences",args.map_sequences, config, default_dict)
    config["map_sequences"] = map_sequences

    ## map_info
    map_info = qcfunk.check_arg_config_default("map_info",args.map_info, config, default_dict)
    config["map_info"] = map_info

    ## input_crs
    input_crs = qcfunk.check_arg_config_default("input_crs",args.input_crs, config, default_dict)
    config["input_crs"] = input_crs

    ## colour_map_by
    colour_map_by = qcfunk.check_arg_config_default("colour_map_by",args.colour_map_by, config, default_dict)
    config["colour_map_by"] = colour_map_by



def report_group_to_config(args,config,default_dict):
    ## sequencing_centre
    sequencing_centre = qcfunk.check_arg_config_default("sequencing_centre",args.sequencing_centre, config, default_dict)
    config["sequencing_centre"] = sequencing_centre

    ## display_name
    display_name = qcfunk.check_arg_config_default("display_name", args.display_name, config, default_dict)
    config["display_name"] = display_name
    
    ## colour_by
    colour_by = qcfunk.check_arg_config_default("colour_by",args.colour_by, config, default_dict)
    config["colour_by"] = colour_by

    ## tree_fields
    tree_fields = qcfunk.check_arg_config_default("tree_fields",args.tree_fields, config, default_dict)
    config["tree_fields"] = tree_fields

    ## label_fields
    label_fields = qcfunk.check_arg_config_default("label_fields",args.label_fields, config, default_dict)
    if not label_fields:
        config["label_fields"] = False

    ##date_fields
    date_fields = qcfunk.check_arg_config_default("date_fields", args.date_fields, config, default_dict)
    config["date_fields"] = date_fields

    ##sample date column
    sample_date_column = qcfunk.check_arg_config_default("sample_date_column", args.sample_date_column,config,default_dict)
    config["sample_date_column"] = sample_date_column

    ## node-summary
    node_summary = qcfunk.check_arg_config_default("node_summary",args.node_summary, config, default_dict)
    config["node_summary"] = node_summary

    ## table_fields
    table_fields = qcfunk.check_arg_config_default("table_fields",args.table_fields, config, default_dict)
    config["table_fields"] = table_fields

    ## include_snp_table
    include_snp_table = qcfunk.check_arg_config_default("include_snp_table",args.include_snp_table, config, default_dict)
    config["include_snp_table"] = include_snp_table

    ## include_bars
    include_bars = qcfunk.check_arg_config_default("include_bars",args.include_bars, config, default_dict)
    config["include_bars"] = include_bars

    ## omit-appendix
    omit_appendix = qcfunk.check_arg_config_default("omit_appendix",args.omit_appendix, config, default_dict)
    config["omit_appendix"] = omit_appendix

    ## no-snipit
    no_snipit = qcfunk.check_arg_config_default("no_snipit",args.no_snipit, config, default_dict)
    config["no_snipit"] = no_snipit
    
    ## private
    private = qcfunk.check_arg_config_default("private",args.private, config, default_dict)
    config["private"] = private

def make_full_civet_table(query_dict, tree_fields, label_fields, input_column, outdir, table_fields, include_snp_table):

    df_dict_incog = defaultdict(list)
    df_dict_seqprovided = defaultdict(list)

    incog = 0
    seqprovided = 0
    incogs = False
    seqprovideds = False

    for query in query_dict.values():

        if query.in_db:
            df_dict = df_dict_incog
            incog += 1
        else:
            df_dict = df_dict_seqprovided
            seqprovided += 1
        
        if query.display_name != query.name:
            df_dict["Query ID"].append(query.display_name.replace("|","\|"))
        else:
            df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        if query.in_db: 
            df_dict["Sequence name in Tree"].append(query.name)        

        df_dict["Sample date"].append(query.sample_date)

        if not query.in_db: 
            df_dict["Closest sequence in Tree"].append(query.closest)
            df_dict["Distance to closest sequence"].append(query.closest_distance)
            df_dict["SNPs"].append(query.snps)

        df_dict["UK lineage"].append(query.uk_lin)
        df_dict["Global lineage"].append(query.global_lin)
        df_dict["Phylotype"].append(query.phylotype)

        if query.tree != "NA":
            tree_number = query.tree.split("_")[-1]
            pretty_tree = "Tree " + str(tree_number)
            df_dict["Tree"].append(pretty_tree)
        else:
            df_dict["Tree"].append("NA") #this should never happen, it's more error catching

        if tree_fields != []:
            for i in tree_fields:
                df_dict[i].append(query.attribute_dict[i])
        
        if label_fields != []:
            for i in label_fields: 
                if i not in tree_fields and i != "sample_date" and i != input_column:
                    df_dict[i].append(query.attribute_dict[i])

    if incog != 0:
        df_incog = pd.DataFrame(df_dict_incog)
        file_name = os.path.join(outdir,"Sequences_already_in_cog")
        df_incog.to_csv(file_name, index=False)
        df_incog.set_index("Query ID", inplace=True)
        incogs = True
    
    if seqprovided != 0:
        df_seqprovided = pd.DataFrame(df_dict_seqprovided)
        file_name = os.path.join(outdir,"Sequences_provided")
        df_seqprovided.to_csv(file_name, index=False)
        df_seqprovided.set_index("Query ID", inplace=True)
        seqprovideds = True

    output = table_func.make_custom_table(query_dict, table_fields, include_snp_table)

    # print(output)
    # print(len(output))

    # for i in output:
    #     print(type(i))


    # for i in output:
    #     print("item")
    #     print(i)

    return output