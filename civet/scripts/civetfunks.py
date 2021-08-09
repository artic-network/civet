#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import date
import datetime as dt
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

def get_defaults():
    default_dict = {
                    "title": "## Cluster investigation",
                    "outbreak_id": "",
                    "report_date": today,# date investigation was opened
                    "authors": "", # List of authors, affiliations and contact details
                    "description": "",
                    "conclusions": "",
                    "max_ambiguity":0.5,
                    "min_length":10000,
                    "num_seqs":0,
                    "no_temp":False,
                    "datadir":"",
                    "input_column":"name",
                    "data_column":"central_sample_id",
                    "database_sample_date_column":"sample_date",
                    "sample_date_column":"sample_date",
                    "display_name":False,
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
                    "from_metadata":False,
                    "date_restriction":False,
                    "date_range_start":False,
                    "date_range_end":False,
                    "launch_browser":False,
                    "background_metadata_all":False,
                    "node_summary":"country",
                    "date_window":7,
                    "colour_by":"adm1=Paired",
                    "label_fields":False,
                    "date_fields":False,
                    "remote":False,
                    "no_snipit":False,
                    "global_snipit": False,
                    "reference_name": "Reference",
                    "remove_snp_table":False,
                    "include_bars":False,
                    "omit_trees":False,
                    "omit_appendix":True,
                    "table_fields":["sample_date", "uk_lineage", "lineage", "phylotype"],
                    "context_table_summary":False,
                    "threads":1,
                    "force":True,
                    "cluster":False,
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674,
                    "protect": False,
                    "output_prefix":"civet",
                    "update":False,
                    "safety_level":1
                    }
    return default_dict

def define_seq_db(config):
    config["seq_db"] = config["background_seqs"]
    

def check_adm2_values(config):

    get_acceptable_adm2(config)
    accepted_adm2 = config["clean_locs"]

    with open(config["query"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if "adm2" in header:
            for row in reader:
                if row["adm2"].upper().replace(" ","_") not in accepted_adm2 and "|" not in row["adm2"] and row["adm2"] != "": 
                    adm2_value = row["adm2"]
                    sys.stderr.write(qcfunk.cyan(f'Error: {adm2_value} not a valid adm2 region.\n Find a list of valid adm2 values at:\nhttps://artic-network.github.io/civet/geographic_data.html\n'))
                    sys.stderr.write(qcfunk.cyan(f'Please note: if you have a region that encompasses multiple adm2 regions eg West Midlands, list the correct adm2 regions separated by the "|" symbol to indicate ambiguity.\n'))
                    sys.exit(-1)
                elif "|" in row["adm2"]:
                    adm2_value = row["adm2"].split("|")
                    for i in adm2_value:
                        if i not in accepted_adm2:
                            sys.stderr.write(qcfunk.cyan(f'Error: {i} found in the ambiguity code {row["adm2"]} not a valid adm2 region.\n Find a list of valid adm2 values at:\nhttps://artic-network.github.io/civet/geographic_data.html\n'))
                            sys.exit(-1)


def get_package_data(thisdir,config):
    reference_fasta = pkg_resources.resource_filename('civet', 'data/reference.fasta')
    outgroup_fasta = pkg_resources.resource_filename('civet', 'data/outgroup.fasta')
    polytomy_figure = pkg_resources.resource_filename('civet', 'data/polytomies.png')
    report_args = pkg_resources.resource_filename('civet', 'data/report_arguments.txt')
    footer_fig = pkg_resources.resource_filename('civet', 'data/footer.png')
    clean_locs_file = pkg_resources.resource_filename('civet', 'data/mapping_files/adm2_cleaning.csv')
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
    
    config["clean_locs_file"] = clean_locs_file
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
    - cog_global_2020-XX-YY_tree.newick\n\
    - cog_global_2020-XX-YY_metadata.csv\n\
    - cog_global_2020-XX-YY_alignment.fasta\n\n"
    +qcfunk.cyan("Please also check that the data directory is correctly specified.\n\n")
    +qcfunk.cyan(f"\
To run civet please either\n1) ssh into CLIMB and run with --CLIMB flag\n\
2) Run using `--remote` flag and your CLIMB username specified e.g. `-uun climb-covid19-smithj`\n\
3) Specify a local directory with the appropriate files, optionally supply a custom metadata file, custom background tree or custom background fasta file\n\n"""))

def rsync_data_from_climb(uun, data_dir):
    rsync_command = f"rsync -avzh --exclude 'cog' --delete-after {uun}@bham.covid19.climb.ac.uk:/cephfs/covid/bham/results/phylogenetics/latest/civet/ '{data_dir}'"
    print(qcfunk.green(f"Syncing civet data to {data_dir}"))
    status = os.system(rsync_command)
    if status != 0:
        sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your user name is a valid CLIMB username e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and are in the UK\n\n"))
        sys.exit(-1)

def get_background_files(data_dir,background_metadata, background_tree, background_seqs,background_metadata_all=False):
    # background_seqs = ""
    # background_tree = ""
    data_date = ""
    for r,d,f in os.walk(data_dir):
        for fn in f:
            if background_seqs == ""  and fn.endswith(".fasta") and fn.startswith("cog_global_"):
                background_seqs = os.path.join(data_dir, fn)
                data_date = fn.split("_")[2]
                if not data_date.startswith("20"):
                    data_date = ""
            elif background_tree == ""  and fn.endswith(".newick") and fn.startswith("cog_global_"):
                background_tree = os.path.join(data_dir, fn)
            elif background_metadata_all and fn.endswith("all.csv") and fn.startswith("cog_global_"):
                background_metadata_all = os.path.join(data_dir, fn)
            elif background_metadata == "" and fn.endswith(".csv") and fn.startswith("cog_global_"):
                background_metadata = os.path.join(data_dir, fn)
                
    return background_seqs, background_tree, background_metadata, data_date, background_metadata_all
    

def get_remote_data(uun,background_metadata,background_trees, background_seqs,data_dir,config):
    config["remote"]= True
    if uun:
        config["username"] = uun
        rsync_data_from_climb(uun, data_dir)
    elif "username" in config:
        uun = config["username"]
        rsync_data_from_climb(uun, data_dir)
    elif "uun" in config:
        uun = config["uun"]
        rsync_data_from_climb(uun, data_dir)
    else:
        rsync_command = f"rsync -avzh --exclude 'cog' --delete-after  bham.covid19.climb.ac.uk:/cephfs/covid/bham/results/phylogenetics/latest/civet/ '{data_dir}'"
        print(f"Syncing civet data to {data_dir}")
        status = os.system(rsync_command)
        if status != 0:
            sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your ssh is configured with Host bham.covid19.climb.ac.uk\nAlternatively enter your CLIMB username with -uun e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and check if you are in the UK\n\n"))
            sys.exit(-1)

    background_seqs, background_tree, background_metadata, data_date, background_metadata_all = get_background_files(data_dir,background_metadata, background_trees, background_seqs)

    config["datadir"] = data_dir
    config["data_date"] = data_date
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

        print(qcfunk.green("Found data:"))
        print("    -",background_seqs)
        print("    -",background_metadata)
        print("    -",background_tree,"\n")

def get_datadir(args_climb,args_uun,args_datadir,args_metadata, args_tree, args_fasta, cwd,config):
    data_dir = ""
    background_metadata = ""
    background_seqs = ""
    background_tree = ""
    remote= config["remote"]
    cog_all = False

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
    if args_tree:
        expanded_path = os.path.expanduser(args_tree)
        background_tree= os.path.join(cwd, expanded_path)
        if not os.path.exists(background_tree):
            sys.stderr.write(qcfunk.cyan(f"Error: can't find tree file at {background_tree}.\n"))
            sys.exit(-1)

    elif "background_tree" in config:
        if config["background_tree"]:
            expanded_path = os.path.expanduser(config["background_tree"])
            background_tree= os.path.join(config["path_to_query"], expanded_path)
            if not os.path.exists(background_tree):
                sys.stderr.write(qcfunk.cyan(f"Error: can't find tree file at {background_tree}.\n"))
                sys.exit(-1)

    if args_fasta:
        expanded_path = os.path.expanduser(args_fasta)
        background_seqs = os.path.join(cwd, expanded_path)
        if not os.path.exists(background_seqs):
            sys.stderr.write(qcfunk.cyan(f"Error: can't find metadata file at {background_seqs}.\n"))
            sys.exit(-1)

    elif "background_sequences" in config:
        if config["background_sequences"]:
            expanded_path = os.path.expanduser(config["background_sequences"])
            background_seqs = os.path.join(config["path_to_query"], expanded_path)
            if not os.path.exists(background_seqs):
                sys.stderr.write(qcfunk.cyan(f"Error: can't find fasta file at {background_seqs}.\n"))
                sys.exit(-1)
            
    if args_climb:
        data_dir = "/cephfs/covid/bham/results/phylogenetics/latest/civet/cog"
        cog_all = False
        if os.path.exists(data_dir):
            config["remote"] = False
            config["username"] = ""
        else:
            sys.stderr.write(qcfunk.cyan(f"Error: --CLIMB argument called, but CLIMB data path doesn't exist.\n"))
            sys.exit(-1)

    elif args_datadir:
        data_dir = os.path.join(cwd, args_datadir)

    elif "datadir" in config:
        if config["datadir"]:
            expanded_path = os.path.expanduser(config["datadir"])
            data_dir = os.path.join(config["path_to_query"], expanded_path)
        else:
            data_dir = os.path.join(cwd, "civet-cat")


    if not remote:
        if not os.path.exists(data_dir):
            print_data_error(data_dir)
            sys.exit(-1)
        
        
        background_seqs, background_tree, background_metadata, data_date, background_metadata_all = get_background_files(data_dir,background_metadata, background_tree,background_seqs, cog_all)
        
        
        config["datadir"] = data_dir
        config["data_date"] = data_date

        if not os.path.isfile(background_tree) or not os.path.isfile(background_seqs) or not os.path.isfile(background_metadata):
            print_data_error(data_dir)
            sys.exit(-1)
        else:
            config["background_metadata"] = background_metadata
            config["background_seqs"] = background_seqs
            config["background_tree"] = background_tree
            config["background_metadata_all"] = background_metadata_all

            print("Found data:")
            print("    -",background_seqs)
            print("    -",background_metadata)
            print("    -",background_metadata_all)
            print("    -",background_tree,"\n")

    elif remote:
        
        get_remote_data(args_uun, background_metadata, background_tree, background_seqs, data_dir, config)

    config["datadir"]=data_dir

def check_update_dependencies(config):
    if "from_metadata" in config:
        if not config["from_metadata"]:
            sys.stderr.write(qcfunk.cyan('Error: `--from-metadata` search term required to run in `update` mode\n'))
            sys.exit(-1)
    else:
        sys.stderr.write(qcfunk.cyan('Error: `--from-metadata` search term required to run in `update` mode\n'))
        sys.exit(-1)

def configure_update(update_arg,udpate_arg,config):
    if udpate_arg:
        config["update"] = True
    qcfunk.add_arg_to_config("update",update_arg, config)
    if config["update"]:
        check_update_dependencies(config)
        for i in ["treedir","tempdir","summary_dir","sequencing_centre_dest","seq_db","rel_outdir","path_to_query","outfile","outdir","background_metadata","background_seqs","background_tree","collapse_summary","filtered_background_metadata"]:
            config[i]=""
        config["colour_by"]="new:Paired"
        config["tree_fields"]="new"

        config["table_fields"]=list(config["table_fields"]).append("new")


def check_cluster_dependencies(config):
    if not "query" in config:
        sys.stderr.write(qcfunk.cyan('Error: input.csv required to run `cluster` civet\n'))
        sys.exit(-1)
    if config["update"]:
        sys.stderr.write(qcfunk.cyan('Error: specify one of either `cluster` or `update`\n'))
        sys.exit(-1)

def configure_cluster(config):
    if config["cluster"]:
        check_cluster_dependencies(config)
    config["colour_by"]="new:Paired"
    config["table_fields"]=list(config["table_fields"]).append("new")
    config["down_distance"]=100


def check_for_update(update_query,config):

    check_update_dependencies(config)
    temp_query = os.path.join(config["outdir"], "pre_update_query.csv")
    new_query = qcfunk.generate_query_from_metadata(temp_query,config["from_metadata"],config["background_metadata"],config)

    old_query = []

    with open(config["query"],"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            old_query.append(row[config["input_column"]])

    update_files = False
    new_seqs = []

    with open(update_query,"w") as fw:
        with open(temp_query,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames

            if "new" not in header:
                header.append("new")

            writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
            writer.writeheader()

            for row in reader:
                new_row = row
                if row[config["input_column"]] not in old_query:
                    update_files = True
                    new_row["new"]="True"
                    writer.writerow(new_row)
                else:
                    new_row["new"]="False"
                    writer.writerow(new_row)
    
    if new_seqs:
        from_metadata = config["from_metadata"]
        print(qcfunk.green(f"New sequences identified with {from_metadata} search"))
 
    config["query"] = update_query
    return update_files

def check_for_new_in_cluster(config):
    new_count = 0
    prefix = config["output_prefix"]
    background_metadata = config["background_metadata"]
    cluster_csv = os.path.join(config["outdir"],f"{prefix}.csv")
    with open(cluster_csv, "r") as f:
        reader = csv.DictReader(f)
        if not "new" in reader.fieldnames:
            sys.stderr.write(qcfunk.cyan('Error: `cluster` civet has not run, require column `new` in csv\n'))
            sys.exit(-1)
        for row in reader:
            if row["new"] == "True":
                new_count +=1
    return new_count, cluster_csv


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

def local_lineages_qc(config):

    query_file = config["query"]
    date_format = "%Y-%m-%d"

    if config["local_lineages"]:

        if "adm2" not in config["background_metadata_header"]:
            sys.stderr.write(qcfunk.cyan('Error: no geographic information found for local lineage analysis. Please provide a column in the background metadata with the header "adm2"\n'))
            sys.exit(-1)
        elif "uk_lineage" not in config["background_metadata_header"]:
            sys.stderr.write(qcfunk.cyan('Error: no uk lineage information found for local lineage analysis. Please provide a column in the background metadata with the header "uk_lineage"\n'))
            sys.exit(-1)

        if config["date_restriction"]:
            if config["date_range_start"] and type(config["date_range_start"]) == str:
                try:
                    check_date = dt.datetime.strptime(config["date_range_start"], date_format).date()
                except:
                    print(qcfunk.cyan(f'date range start in incorrect format. Please use i.e. YYYY-MM-DD'))
                    sys.exit(-1)
                
            if config["date_range_end"] and type(config["date_range_end"]) == str:
                try:
                    check_date = dt.datetime.strptime(config["date_range_end"], date_format).date()
                except:
                    print(qcfunk.cyan(f'date range end in incorrect format. Please use i.e. YYYY-MM-DD'))
                    sys.exit(-1)

            if config["date_range_start"] and config["date_range_end"]:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_range_start']} to {config['date_range_end']}"))
            elif config["date_range_start"]:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_range_start']} to present"))
            else:
                print(qcfunk.green(f"Local lineage analysis restricted to {config['date_window']} days around the sampling range"))

        elif config['date_range_start'] or config["date_range_end"]:
            print(qcfunk.cyan("Date restriction data provided but --date-restriction flag not used. Please use --date-restriction flag in config or command line."))
            sys.exit(-1)

        else:
            print(qcfunk.green(f"Local lineage analysis not restricted by time, will show background lineage composition for the whole of the epidemic"))

        

def local_lineages_to_config(central, neighbouring, region, config):

    if config["local_lineages"] == True:
        lineage_tables = []
        for r,d,f in os.walk(os.path.join(config["outdir"],"report", 'figures')):
            for fn in f:
                if fn.endswith("_lineageTable.md"):
                    lineage_tables.append(os.path.join(config["outdir"],"report", 'figures', fn))

        config["lineage_tables"] = lineage_tables
        config["lineage_maps"] = [central, neighbouring, region]
    else:
        config["lineage_tables"] = []
        config["lineage_maps"] = []

def map_sequences_config(config):
    
    background_headers = config["background_metadata_header"]
    query_headers = config["query_metadata_header"]

    if config["map_sequences"]:

        map_inputs = ""
        if config["map_info"]:
            map_inputs = config["map_info"].replace(" ","")
        else:
            sys.stderr.write(qcfunk.cyan('Error: coordinates or outer postcode or adm2 not supplied for mapping sequences. Please provide either x and y columns as a comma separated string, or column header containing outer postcode.'))
            sys.exit(-1)

        if len(map_inputs.split(",")) == 2: #If x and y coordinates are provided
            if not config["input_crs"]:
                sys.stderr.write(qcfunk.cyan('Error: input coordinate system not provided for mapping. Please provide --input-crs eg EPSG:3395'))
                sys.exit(-1)
        else: #If an outer postcode column is provided
            config["input_crs"] = "EPSG:4326"

        relevant_cols = map_inputs.split(",")

        if config["colour_map_by"]:
            relevant_cols.append(config["colour_map_by"])
        
        for map_arg in relevant_cols:
            map_arg = map_arg.replace(" ","")
            if map_arg not in query_headers and map_arg not in background_headers:
                sys.stderr.write(qcfunk.cyan(f"Error: {map_arg} column not found in metadata file or background database for mapping sequences"))
                sys.exit(-1)

        if config["colour_map_by"]:
            if map_inputs == "adm2":
                print(qcfunk.cyan(f"NOTE: --colour-map-by not set up to colour by adm2. Please provide outer postcode or coordinates"))
            else:
                print(qcfunk.green(f"Colouring map by: {config['colour_map_by']}"))
                        

def get_sequencing_centre_header(config):
    
    sc_list = ["PHEC", 'LIVE', 'BIRM', 'PHWC', 'CAMB', 'NORW', 'GLAS', 'EDIN', 'SHEF',
                'EXET', 'NOTT', 'PORT', 'OXON', 'NORT', 'NIRE', 'GSTT', 'LOND', 'SANG',"NIRE"]

    sequencing_centre = config["sequencing_centre"]
    if sequencing_centre in sc_list or sequencing_centre == "DEFAULT":
        package_png = os.path.join("data","headers",f"{sequencing_centre}.png")
        sequencing_centre_source = pkg_resources.resource_filename('civet', package_png)
        print(qcfunk.green(f"Using header file from:") + f" {package_png}\n")
        config["sequencing_centre_source"] = sequencing_centre_source
        config["sequencing_centre_dest"] = os.path.join(config["outdir"],"report","figures",f"{sequencing_centre}.png")
        config["sequencing_centre_file"] = os.path.join(".","figures",f"{sequencing_centre}.png")
        config["sequencing_centre"] = sequencing_centre
    else:
        sc_string = "\n".join(sc_list)
        sys.stderr.write(qcfunk.cyan(f'Error: sequencing centre must be one of the following:\n{sc_string}\n'))
        sys.exit(-1)

def map_group_to_config(args,config):

    ## local_lineages
    qcfunk.add_arg_to_config("local_lineages",args.local_lineages, config)

    ## date_restriction
    qcfunk.add_arg_to_config("date_restriction",args.date_restriction, config)
    ## date_range_start
    qcfunk.add_arg_to_config("date_range_start",args.date_range_start, config)

    ## date_range_end
    qcfunk.add_arg_to_config("date_range_end",args.date_range_end, config)

    ## date_window
    qcfunk.add_arg_to_config("date_window",args.date_window, config)

    ## map_sequences
    qcfunk.add_arg_to_config("map_sequences",args.map_sequences, config)

    ## map_info
    qcfunk.add_arg_to_config("map_info",args.map_info, config)
    ## input_crs
    qcfunk.add_arg_to_config("input_crs",args.input_crs, config)

    ## colour_map_by
    qcfunk.add_arg_to_config("colour_map_by",args.colour_map_by, config)


def report_group_to_config(args,config):
    ## sequencing_centre
    qcfunk.add_arg_to_config("sequencing_centre",args.sequencing_centre, config)

    ## display_name
    qcfunk.add_arg_to_config("display_name", args.display_name, config)
    
    ## colour_by
    qcfunk.add_arg_to_config("colour_by",args.colour_by, config)

    ## tree_fields
    qcfunk.add_arg_to_config("tree_fields",args.tree_fields, config)

    ## label_fields
    qcfunk.add_arg_to_config("label_fields",args.label_fields, config)

    ##date_fields
    qcfunk.add_arg_to_config("date_fields", args.date_fields, config)

    ##sample date column
    qcfunk.add_arg_to_config("sample_date_column", args.sample_date_column,config)
    qcfunk.add_arg_to_config("database_sample_date_column", args.database_sample_date_column, config)

    ## node-summary
    qcfunk.add_arg_to_config("node_summary",args.node_summary, config)

    ## table_fields
    qcfunk.add_arg_to_config("table_fields",args.table_fields, config)
    ## remove_snp_table
    qcfunk.add_arg_to_config("remove_snp_table",args.remove_snp_table, config)

    ## include_bars
    qcfunk.add_arg_to_config("include_bars",args.include_bars, config)

    ## omit-appendix
    qcfunk.add_arg_to_config("omit_appendix",args.omit_appendix, config)

    ## no-snipit
    qcfunk.add_arg_to_config("no_snipit",args.no_snipit, config)

    ## global snipit figures

    qcfunk.add_arg_to_config("global_snipit", args.global_snipit, config)

    ## omit-trees
    qcfunk.add_arg_to_config("omit_trees", args.omit_trees, config)

    ##context-table
    qcfunk.add_arg_to_config("context_table_summary", args.context_table_summary, config)
    

def make_full_civet_table(query_dict, full_taxon_dict, tree_fields, label_fields, input_column, outdir, table_fields):

    df_dict = defaultdict(list)

    for name,taxon in full_taxon_dict.items():
            
        if name in query_dict or taxon.protected:

            df_dict["Query ID"].append(taxon.name.replace("|","\|"))
            df_dict["Name used in report"].append(taxon.display_name.replace("|","\|"))

            if taxon.in_db: 
                df_dict["sequence_name"].append(taxon.name)   
            else:
                df_dict["sequence_name"].append("")     

            df_dict["sample_date"].append(taxon.sample_date)

            if not taxon.in_db and not taxon.protected: 
                df_dict["closest_sequence"].append(taxon.closest)
                df_dict["distance_to_closest"].append(taxon.closest_distance)
                df_dict["SNPs"].append(taxon.snps)
            else:
                df_dict["closest_sequence"].append("")
                df_dict["distance_to_closest"].append("")
                df_dict["SNPs"].append("")

            if taxon.in_db:
                df_dict["in_cog"].append("True")
            elif taxon.protected:
                df_dict["in_cog"].append("Background sequence")
            else:
                df_dict["in_cog"].append("False")

            # df_dict["UK_lineage"].append(taxon.uk_lineage)
            # df_dict["lineage"].append(taxon.global_lineage)
            # df_dict["phylotype"].append(taxon.phylotype)

            if taxon.tree != "NA":
                tree_number = taxon.tree.split("_")[-1]
                pretty_tree = "Tree " + str(tree_number)
                df_dict["Tree"].append(pretty_tree)
            else:
                df_dict["Tree"].append("") #this should never happen, it's more error catching

            if tree_fields != []:
                for i in tree_fields:
                    df_dict[i].append(taxon.attribute_dict[i])
            
            if label_fields != []:
                for i in label_fields: 
                    if i not in tree_fields and i != "sample_date" and i != input_column:
                        df_dict[i].append(taxon.attribute_dict[i])

    df = pd.DataFrame(df_dict)

    file_name = os.path.join(outdir,"civet_metadata.csv")
    df.to_csv(file_name, index=False)


def anonymise_sequences(taxon_dict, query_dict, safety_level, from_metadata): 
    #if it's in the query and the display_name is given, then that should be what the seq name is

    count = 0
    for name,tax in sorted(taxon_dict.items(), key=lambda x: random.random()):
        
        if (name in query_dict and not from_metadata): 
            tax.display_name = tax.input_display_name #ie don't anonymise it if they've provided it themselves OR safe status is 0

        else:
            if safety_level == "0" or safety_level == "2":
                tax.display_name = tax.input_display_name

            elif safety_level == "1":
                if (name not in query_dict or from_metadata != "") and tax.country  == "UK":
                    display_name = "seq_" + str(count)
                    count += 1
                    tax.display_name = display_name

                elif tax.country != "UK":
                    tax.display_name = tax.input_display_name

        taxon_dict[name] = tax

    return taxon_dict


def generate_labels(tax,safety_level, custom_tip_fields):
    
    name = tax.display_name
    date = tax.sample_date
    
    display_name = f"{name}|{date}"
    
    if "location_label" in tax.attribute_dict.keys() and safety_level != "2": #if it's being run locally OR if safe status is on no adm2 for distribution
        adm2 = tax.attribute_dict["location_label"]
        display_name = f"{name}|{adm2}|{date}"

    count = 0
    if len(custom_tip_fields) > 0: 
        for label_element in custom_tip_fields:
            if label_element == "adm2":
                label_element = tax.attribute_dict["location_label"]
            else:   
                display_name = display_name + "|" + tax.attribute_dict[label_element]
            count += 1
    
    return display_name

def get_acceptable_adm2(config):

    GADM_adm2 = [
    'BARNSLEY', 'BATH_AND_NORTH_EAST_SOMERSET', 'BEDFORDSHIRE', 'BIRMINGHAM', 'BLACKBURN_WITH_DARWEN', 'BLACKPOOL', 'BOLTON', 'BOURNEMOUTH', 'BRACKNELL_FOREST', 'BRADFORD', 'BRIGHTON_AND_HOVE', 'BRISTOL', 'BUCKINGHAMSHIRE', 'BURY',
    'CALDERDALE', 'CAMBRIDGESHIRE', 'CENTRAL_BEDFORDSHIRE', 'CHESHIRE_EAST', 'CHESHIRE_WEST_AND_CHESTER', 'CORNWALL', 'COVENTRY', 'CUMBRIA', 
    'DARLINGTON', 'DERBY', 'DERBYSHIRE', 'DEVON', 'DONCASTER', 'DORSET', 'DUDLEY', 'DURHAM', 
    'EAST_RIDING_OF_YORKSHIRE', 'EAST_SUSSEX', 'ESSEX', 
    'GATESHEAD', 'GLOUCESTERSHIRE', 'GREATER_LONDON', 
    'HALTON', 'HAMPSHIRE', 'HARTLEPOOL', 'HEREFORDSHIRE', 'HERTFORDSHIRE', 
    'ISLE_OF_WIGHT', 'ISLES_OF_SCILLY', 
    'KENT', 'KINGSTON_UPON_HULL', 'KIRKLEES', 'KNOWSLEY', 
    'LANCASHIRE', 'LEEDS', 'LEICESTER', 'LEICESTERSHIRE', 'LINCOLNSHIRE', 'LUTON', 
    'MANCHESTER', 'MEDWAY', 'MIDDLESBROUGH', 'MILTON_KEYNES', 
    'NEWCASTLE_UPON_TYNE', 'NORFOLK', 'NORTH_LINCOLNSHIRE', 'NORTH_SOMERSET', 'NORTH_TYNESIDE', 'NORTH_YORKSHIRE', 'NORTHAMPTONSHIRE', 'NORTHUMBERLAND', 'NOTTINGHAM', 'NOTTINGHAMSHIRE', 
    'OLDHAM', 'OXFORDSHIRE', 
    'PETERBOROUGH', 'PLYMOUTH', 'POOLE', 'PORTSMOUTH', 
    'READING', 'REDCAR_AND_CLEVELAND', 'ROCHDALE', 'ROTHERHAM', 'RUTLAND', 
    'SAINT_HELENS', 'SALFORD', 'SANDWELL', 'SEFTON', 'SHEFFIELD', 'SHROPSHIRE', 'SLOUGH', 'SOLIHULL', 'SOMERSET', 'SOUTH_GLOUCESTERSHIRE', 'SOUTH_TYNESIDE', 'SOUTHAMPTON', 'SOUTHEND-ON-SEA', 'STAFFORDSHIRE', 'STOCKPORT', 'STOCKTON-ON-TEES', 'STOKE-ON-TRENT', 'SUFFOLK', 'SUNDERLAND', 'SURREY', 'SWINDON', 
    'TAMESIDE', 'TELFORD_AND_WREKIN', 'THURROCK', 'TORBAY', 'TRAFFORD', 'WAKEFIELD', 'WALSALL', 'WARRINGTON', 'WARWICKSHIRE', 'WEST_BERKSHIRE', 'WEST_SUSSEX', 'WIGAN', 'WILTSHIRE', 'WINDSOR_AND_MAIDENHEAD', 'WIRRAL', 'WOKINGHAM', 'WOLVERHAMPTON', 'WORCESTERSHIRE', 'YORK',
    'ANTRIM_AND_NEWTOWNABBEY', 'ARMAGH_BANBRIDGE_AND_CRAIGAVON', 'BELFAST', 'CAUSEWAY_COAST_AND_GLENS', 'DERRY_AND_STRABANE', 'FERMANAGH_AND_OMAGH', 'LISBURN_AND_CASTLEREAGH', 'MID_AND_EAST_ANTRIM', 'MID_ULSTER', 'NEWRY_MOURNE_AND_DOWN', 'NORTH_DOWN_AND_ARDS', 'TYRONE', 'ANTRIM', 'ARMAGH', 'FERMANAGH', 'LONDONDERRY', 'DOWN',
    'ABERDEEN', 'ABERDEENSHIRE', 'ANGUS', 'ARGYLL_AND_BUTE', 'CLACKMANNANSHIRE', 'DUMFRIES_AND_GALLOWAY', 'DUNDEE', 'EAST_AYRSHIRE', 'EAST_DUNBARTONSHIRE', 'EAST_LOTHIAN', 'EAST_RENFREWSHIRE', 'EDINBURGH', 'EILEAN_SIAR', 'FALKIRK', 'FIFE', 
    'GLASGOW', 'HIGHLAND', 'INVERCLYDE', 'MIDLOTHIAN', 'MORAY', 'NORTH_AYRSHIRE', 'NORTH_LANARKSHIRE', 'ORKNEY_ISLANDS', 'PERTHSHIRE_AND_KINROSS', 'RENFREWSHIRE', 'SCOTTISH_BORDERS', 'SHETLAND_ISLANDS', 'SOUTH_AYRSHIRE', 'SOUTH_LANARKSHIRE', 'STIRLING', 'WEST_DUNBARTONSHIRE', 'WEST_LOTHIAN',
    'ANGLESEY', 'BLAENAU_GWENT', 'BRIDGEND', 'CAERPHILLY', 'CARDIFF', 'CARMARTHENSHIRE', 'CEREDIGION', 'CONWY', 'DENBIGHSHIRE', 'FLINTSHIRE', 'GWYNEDD', 'MERTHYR_TYDFIL', 'MONMOUTHSHIRE', 'NEATH_PORT_TALBOT', 'NEWPORT', 'PEMBROKESHIRE', 'POWYS', 'RHONDDA_CYNON_TAFF', 'SWANSEA', 'TORFAEN', 'VALE_OF_GLAMORGAN', 'WREXHAM',
    'GUERNSEY', "JERSEY"
    ]

    config["clean_locs"] = GADM_adm2



def header(v):
    print(qcfunk.green("""\n
                                    __              __    
                              ____ |__|__  __ _____/  |_ 
                             / ___\|  \  \/ // __ \   __|
                            \  \___|  |\   /\  ___/|  |  
                             \____/ __| \_/  \____/ __|  

                **** Cluster Investigation & Virus Epidemiology Tool ****
                """)+qcfunk.green(f"""
                                        {v}""")+qcfunk.green("""
                        ****************************************
                                                                
                            Aine O'Toole & Verity Hill       
                                    Rambaut Group              
                                Edinburgh University          
\n"""))

def preamble(v):
    header(v)
    funding()
    acknowledgements()

def funding():
    print(qcfunk.green("""
                    Funding:                
                                                                
                                    ARTIC Network               
                        Wellcome Trust Collaborators Award      
                                    206298/Z/17/Z               
                                                                
                            COVID-19 Genomics UK Consortium     
                        UK Department of Health and Social Care 
                            UK Research and Innovation          
                                                                
                                    ReservoirDOCs               
                    European Research Council Consolidator Grant
                                    ERC-2016-COG                
                                                             
"""))

def acknowledgements():
    print(qcfunk.green("""
                    Code contributors:           
                                                            
                        Ben Jackson         gofasta       
                        JT McCrone          clusterfunk     
                        Stefan Rooke        local map 
                        Andrew Rambaut      jclusterfunk    
                                                            
                    Acknowledgements:            
                                                            
                    We thank the following for helpful suggestions, 
                    comments, beta-testing, feature requests and
                    patience.                
                                                            
                        :nickloman:         :mattloose:     
                        :mattbashton:       :tomconnor:     
                        :rebeccadewar:      :martinmchugh:    
                        :richardmyers:      :meerachand:    
                        :samnicholls:       :radpoplawski:   
                        :davidaanensen:     :benlindsey:    
                        :jeffbarrett:       :derekfairley:   
                        :josephhughes:      :davidrobertson:  
                        :richardorton:      :mattholden:
                        :ulfschaefer:       :nataliegroves:   
                        :nikosmanesis:      :jaynaraghwani:   
"""))

def be_arty():
    logo()

def logo():
    print("""
                                       &@                                       
                           *@@,,,,,,,,,,,,,,,,,,,,,@@/                          
                      %@,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@%                     
                   @,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@                  
                @,,,,,,,,,,,,,,,,%(**(,,,,,,,,,,,,,,,,,,,,,,,,,,@               
              @,,,,,,,,,,,,,,,%(((((%**%((((((%,,,,,,,,,,,,,,,,,,,@             
            @,,,,,,,,,,,,,,,%((((((((((%**%((((%,,,,,,,,,,,,,,,,,,,,@           
          @,,,,,,,,,,,,,,,%,,,,,,%(((######((((%,,,,,,,,,,,,,,,,,,,,,,@         
         @,,,,,,,,,%(((%(****#,,,,###****##((((#,,,,,,,,,,,,,,,,,,,,,,,@        
        @,,,,,,,,,,,,,,,((((%(((((##***%#%,,,,,#,,,,,,,,,,,,,,,,,,,,,,,,@       
       @,,,,,,,,,,,,,,#(((((((((((%((*%(((((#,,*,,,,,,,,,,,,,,,,,,,,,,,,,@      
      /*,,,,,,,,,,,,%((((((((((((((((##%(*********%,,,,,,,,,,,,,,,,,,,,,,*      
      @,,,,,,,,,,,%((((((((((((((((((((####  #(%*****%,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,,,,,,(((((((((((((%##%%########(**%**,,,,,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,/**#((((%%***#%#(%***********##,/##%/**%,,,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,****(((((((((((%*****(, ,,,.*/***%(*/###%,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,%****((((((((((((((((#*****,,,,,,,,,,,*%%,%###%%%/,**,,,,,,,,,@     
       @,,,,%******((((((((((((((((((%***%,,,,,,,,,,###,,,,,,,/#**%,,,,,,@      
       (*,,/*********((((((((((((((((%*,,,,,,,,,,,###,,,,/*,,,*(****,,,,/.      
        (%****(%#****/(((((((((((((,,,,,,,,,,,,,,##(,,,,/(,,,%,,****,,,*#       
          @*******##((((((((((((((((((((*,,,,,,,%#(,,,,,,,,,/,,***%,(%@         
           @***************************%,,,,,,,,%#(,,,,,,,%,,%*,,,%,,@          
             @********************%/*******,,,,,,,%,,,*#/(,*,,,,,,,@            
               @********************/******%,,,,%##%,,,,,,,,,,,,,@              
                 @@*******************###((((%#####,,,,,,,,,,,@@                
                    *@************%##((***((((####%,,,,,,,,@*                   
                         @@(***/#%*****%(((((%###,,,,,@@                        
                               @@@@/*((((((((%@@@@                              

""")
