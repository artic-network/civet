#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import datetime 
import tempfile
import pkg_resources
import yaml

from reportfunk.funks import io_functions as qcfunk

def get_defaults():
    default_dict = {"maxambig":0.5,
                    "minlen":10000,
                    "datadir":"civet-cat",
                    "input_column":"name",
                    "data_column":"central_sample_id",
                    "display_name":None,
                    "private":False,
                    "distance":2,
                    "up_distance":2,
                    "down_distance":2,
                    "sequencing_centre":"DEFAULT",
                    "threshold":1,
                    "add_bars":False,
                    "tree_fields":"adm1",
                    "date_range_start":False,
                    "date_range_end":False,
                    "date_window":7,
                    "global_search":False,
                    "label_fields":"sequence_name",#this was none
                    "date_fields":"sample_date",#this was none
                    "graphic_dict":"adm1",
                    "date_restriction":False,
                    "local_lineages":False,
                    "map_sequences":False,
                    "delay_collapse":False,
                    "table_fields":["sample_date", "uk_lineage", "lineage", "phylotype"],
                    "threads":1,
                    "force":True,
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674
                    }
    return default_dict

def define_seq_db(config,default_dict):
    config["seq_db"] = config["background_seqs"]
    

def get_package_data(cog_report,thisdir,config,default_dict):
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

    if cog_report:
        report_template = os.path.join(thisdir, 'scripts','COG_template.pmd')
    elif "cog_report" in config:
        report_template = os.path.join(thisdir, 'scripts','COG_template.pmd')
    else:
        report_template = os.path.join(thisdir, 'scripts','civet_template.pmd')
    
    if not os.path.exists(report_template):
        sys.stderr.write(qcfunk.cyan(f'Error: cannot find report_template at {report_template}\n'))
        sys.exit(-1)

    config["report_template"] = report_template

def print_data_error():
    sys.stderr.write(qcfunk.cyan(f"Error: data directory not found at {data_dir}.\n")+ f"""The directory should contain the following files:\n\
    - cog_global_tree.nexus\n\
    - cog_global_metadata.csv\n\
    - cog_global_alignment.fasta\n\
To run civet please either\n1) ssh into CLIMB and run with --CLIMB flag\n\
2) Run using `--remote-sync` flag and your CLIMB username specified e.g. `-uun climb-covid19-otoolexyz`\n\
3) Specify a local directory with the appropriate files\n\n""")

def rsync_data_from_climb(uun, data_dir):
    rsync_command = f"rsync -avzh {uun}@bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{data_dir}'"
    print(qcfunk.green(f"Syncing civet data to {data_dir}"))
    status = os.system(rsync_command)
    if status != 0:
        sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your user name is a valid CLIMB username e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and are in the UK\n\n"))
        sys.exit(-1)

def get_remote_data(uun,data_dir,config):
    config["remote"]= True
    if uun:
        config["username"] = uun
        rsync_data_from_climb(uun, data_dir)
    elif "username" in config:
        uun = config["username"]
        rsync_data_from_climb(uun, data_dir)
    else:
        rsync_command = f"rsync -avzh bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{data_dir}'"
        print(f"Syncing civet data to {data_dir}")
        status = os.system(rsync_command)
        if status != 0:
            sys.stderr.write(qcfunk.cyan("Error: rsync command failed.\nCheck your ssh is configured with Host bham.covid19.climb.ac.uk\nAlternatively enter your CLIMB username with -uun e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and check if you are in the UK\n\n"))
            sys.exit(-1)

    background_metadata = ""
    background_seqs = ""
    background_tree = ""

    background_metadata = os.path.join(data_dir,"civet-cat","background_metadata.csv")
    background_seqs= os.path.join(data_dir,"civet-cat","cog_global_alignment.fasta")
    background_tree = os.path.join(data_dir,"civet-cat","cog_global_tree.nexus")

    if not os.path.isfile(background_tree) or not os.path.isfile(background_seqs) or not os.path.isfile(background_metadata):
        print_data_error()
        sys.exit(-1)
    else:
        config["background_metadata"] = background_metadata
        config["background_seqs"] = background_seqs
        config["background_tree"] = background_tree

        print("Found data:")
        print("    -",background_seqs)
        print("    -",background_metadata)
        print("    -",background_tree,"\n")

def get_datadir(args_climb,args_uun,args_datadir,remote,cwd,config,default_dict):
    data_dir = ""
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
        data_dir = os.path.join(cwd, args_datadir)

    else:
        if remote:
            data_dir = cwd
        else:
            data_dir = os.path.join(cwd, default_dict["datadir"])

    if not remote:
        if not os.path.exists(data_dir):
            print_data_error()
            sys.exit(-1)
            
        background_metadata = ""
        background_seqs = ""
        background_tree = ""
        
        background_metadata = os.path.join(data_dir,"cog_global_metadata.csv")
        background_seqs= os.path.join(data_dir,"cog_global_alignment.fasta")
        background_tree = os.path.join(data_dir,"cog_global_tree.nexus")

        if not os.path.isfile(background_tree) or not os.path.isfile(background_seqs) or not os.path.isfile(background_metadata):
            print_data_error()
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
        
        get_remote_data(args_uun, data_dir, config)

    config["datadir"]=data_dir

def local_lineages_qc(config,default_dict):

    query_file = config["query"]

    if config["local_lineages"]:

        if config["date_restriction"]:

            if config["date_range_start"] and type(config["date_range_start"]) == str:
                check_date_format(config["date_range_start"])
            if config["date_range_end"] and type(config["date_range_end"]) == str:
                check_date_format(config["date_range_end"])

            if config["date_range_start"] and config["date_range_end"]:
                print(green(f"Local lineage analysis restricted to {config['date_range_start']} to {config['date_range_end']}"))
            elif config["date_range_start"]:
                print(green(f"Local lineage analysis restricted to {config['date_range_start']} to present"))
            else:
                print(green(f"Local lineage analysis restricted to {config['date_window']} days around the sampling range"))

def get_full_headers():
    #could open and parse the cog metadata here to make it flexible
    headers = ["central_sample_id", "biosample_source_id","sequence_name","secondary_identifier",
                "sample_date","epi_week","country","adm1","adm2","outer_postcode","is_surveillance","is_community",
                "is_hcw","is_travel_history","travel_history","lineage","lineage_support","uk_lineage","acc_lineage","del_lineage","phylotype"]

    return headers

def map_sequences_config(map_sequences,mapping_trait,map_inputs,input_crs,config,default_dict):

    map_settings = False
    query_file = config["query"]
    full_metadata_headers = get_full_headers()

    map_settings = config["map_sequences"]
    
    if map_settings:
        if "map_info" in config:
            map_inputs = config["map_info"].replace(" ","")

        if "mapping_trait" in config:
            mapping_trait = config["mapping_trait"]

        if not map_inputs:
            sys.stderr.write(cyan('Error: coordinates or outer postcode not supplied for mapping sequences. Please provide either x and y columns as a comma separated string, or column header containing outer postcode.'))
            sys.exit(-1)
        else:
            if len(map_inputs.split(",")) == 2: #If x and y coordinates are provided
                if input_crs:
                    crs = input_crs
                elif "input_crs" in config:
                    crs = config["input_crs"]
                else:
                    sys.stderr.write('Error: input coordinate system not provided for mapping. Please provide --input-crs eg EPSG:3395')
                    sys.exit(-1)
            else: #If an outer postcode column is provided        
                crs = "EPSG:4326"
                            
            config["map_info"] = map_inputs.replace(" ","")
            config["input_crs"] = crs

        with open(query_file, newline="") as f:
            reader = csv.DictReader(f)
            column_names = reader.fieldnames
            relevant_cols = []
            map_inputs_lst = map_inputs.split(",")
            for i in map_inputs_lst:
                relevant_cols.append(i)
            
            if mapping_trait:
                relevant_cols.append(mapping_trait)
            
            for map_arg in relevant_cols:
                map_arg = map_arg.replace(" ","")
                if map_arg not in column_names and map_arg not in full_metadata_headers:
                    sys.stderr.write(cyan(f"Error: {map_arg} field not found in metadata file or background database for mapping sequences"))
                    sys.exit(-1)

        if mapping_trait:
            if map_inputs == "adm2":
                print(cyan(f"NOTE: mapping trait provided, but summary map is not designed for showing trait. Please provide more detailed mapping information, eg outer postcode or coordinates"))
            else:
                print(green(f"Colouring map by: " + f"{mapping_trait}"))
            config["mapping_trait"] = mapping_trait
            
        else:
            config["mapping_trait"] = False
            
    else:
        config["map_sequences"] = False
        config["map_info"] = False
        config["input_crs"] = False
        config["mapping_trait"] = False


def get_sequencing_centre_header(config):
    
    sc_list = ["PHEC", 'LIVE', 'BIRM', 'PHWC', 'CAMB', 'NORW', 'GLAS', 'EDIN', 'SHEF',
                'EXET', 'NOTT', 'PORT', 'OXON', 'NORT', 'NIRE', 'GSTT', 'LOND', 'SANG',"NIRE"]


    if sequencing_centre in sc_list or sequencing_centre == "DEFAULT":
        package_png = os.path.join("data","headers",f"{sequencing_centre}.png")
        sequencing_centre_source = pkg_resources.resource_filename('civet', package_png)
        print(green(f"Using header file from:") + f" {package_png}\n")
        config["sequencing_centre_source"] = sequencing_centre_source
        config["sequencing_centre_dest"] = os.path.join(config["outdir"],"figures",f"{sequencing_centre}.png")
        config["sequencing_centre_file"] = os.path.join(".","figures",f"{sequencing_centre}.png")
        config["sequencing_centre"] = sequencing_centre
    else:
        sc_string = "\n".join(sc_list)
        sys.stderr.write(cyan(f'Error: sequencing centre must be one of the following:\n{sc_string}\n'))
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

    ## map_cols
    map_cols = qcfunk.check_arg_config_default("map_cols",args.map_cols, config, default_dict)
    config["map_cols"] = map_cols

    ## input_crs
    input_crs = qcfunk.check_arg_config_default("input_crs",args.input_crs, config, default_dict)
    config["input_crs"] = input_crs

    ## mapping_trait
    mapping_trait = qcfunk.check_arg_config_default("mapping_trait",args.mapping_trait, config, default_dict)
    config["mapping_trait"] = mapping_trait



def report_columns_to_config(args,config,default_dict):
    ## sequencing_centre
    sequencing_centre = qcfunk.check_arg_config_default("sequencing_centre",args.sequencing_centre, config, default_dict)
    config["sequencing_centre"] = sequencing_centre
    
    ## display
    display = qcfunk.check_arg_config_default("display",args.display, config, default_dict)
    config["display"] = display

    ## fields
    fields = qcfunk.check_arg_config_default("fields",args.fields, config, default_dict)
    config["fields"] = fields

    ## label_fields
    label_fields = qcfunk.check_arg_config_default("label_fields",args.label_fields, config, default_dict)
    config["label_fields"] = label_fields

    ## node-summary
    node_summary = qcfunk.check_arg_config_default("node_summary",args.node_summary, config, default_dict)
    config["node_summary"] = node_summary

    ## table_fields
    table_fields = qcfunk.check_arg_config_default("table_fields",args.table_fields, config, default_dict)
    config["table_fields"] = table_fields

    ## snp-data-table
    snp_data_table = qcfunk.check_arg_config_default("snp_data_table",args.snp_data_table, config, default_dict)
    config["snp_data_table"] = snp_data_table

    ## add-bars
    add_bars = qcfunk.check_arg_config_default("add_bars",args.add_bars, config, default_dict)
    config["add_bars"] = add_bars

    ## cog_report
    cog_report = qcfunk.check_arg_config_default("cog_report",args.cog_report, config, default_dict)
    config["cog_report"] = cog_report

    ## omit-appendix
    omit_appendix = qcfunk.check_arg_config_default("omit_appendix",args.omit_appendix, config, default_dict)
    config["omit_appendix"] = omit_appendix
    
    ## private
    private = qcfunk.check_arg_config_default("private",args.private, config, default_dict)
    config["private"] = private

