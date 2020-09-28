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

def define_seq_db(global_arg,config,default_dict):
    global_search = qcfunk.check_arg_config_default("global_search",global_arg, config, default_dict)
    if global_search:
        config["seq_db"] = config["cog_global_seqs"]
    else:
        config["seq_db"] = config["cog_seqs"]


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
    - cog_metadata.csv\n\
    - cog_global_metadata.csv\n\
    - cog_global_alignment.fasta\n\
    - cog_alignment.fasta\n\n\
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

    cog_metadata,cog_global_metadata = ("","")
    cog_seqs = ""
    cog_tree = ""

    cog_seqs = os.path.join(data_dir,"civet-cat","cog_alignment.fasta")
    cog_metadata = os.path.join(data_dir,"civet-cat","cog_metadata.csv")

    cog_global_metadata = os.path.join(data_dir,"civet-cat","cog_global_metadata.csv")
    cog_global_seqs= os.path.join(data_dir,"civet-cat","cog_global_alignment.fasta")

    cog_tree = os.path.join(data_dir,"civet-cat","cog_global_tree.nexus")

    if not os.path.isfile(cog_seqs) or not os.path.isfile(cog_global_seqs) or not os.path.isfile(cog_metadata) or not os.path.isfile(cog_global_metadata) or not os.path.isfile(cog_tree):
        print_data_error()
        sys.exit(-1)

    config["cog_seqs"] = cog_seqs

    config["cog_metadata"] = cog_metadata
    config["cog_global_metadata"] = cog_global_metadata

    config["cog_global_seqs"] = cog_global_seqs
    config["cog_tree"] = cog_tree

    print("Found cog data:")
    print("    -",cog_seqs)
    print("    -",cog_metadata)
    print("    -",cog_global_metadata)
    print("    -",cog_tree,"\n")

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
            
        cog_metadata,cog_global_metadata = ("","")
        cog_seqs = ""
        cog_tree = ""
        
        cog_seqs = os.path.join(data_dir,"cog_alignment.fasta")
        cog_metadata = os.path.join(data_dir,"cog_metadata.csv")

        cog_global_metadata = os.path.join(data_dir,"cog_global_metadata.csv")
        cog_global_seqs= os.path.join(data_dir,"cog_global_alignment.fasta")

        cog_tree = os.path.join(data_dir,"cog_global_tree.nexus")

        if not os.path.isfile(cog_seqs) or not os.path.isfile(cog_global_seqs) or not os.path.isfile(cog_metadata) or not os.path.isfile(cog_global_metadata) or not os.path.isfile(cog_tree):
            print_data_error()
            sys.exit(-1)
        else:
            config["cog_seqs"] = cog_seqs
            config["cog_metadata"] = cog_metadata

            config["cog_global_metadata"] = cog_global_metadata
            config["cog_global_seqs"] = cog_global_seqs
            config["cog_tree"] = cog_tree

            print("Found cog data:")
            print("    -",cog_seqs)
            print("    -",cog_metadata)
            print("    -",cog_global_metadata)
            print("    -",cog_tree,"\n")

    elif remote:
        
        get_remote_data(args_uun, data_dir, config)

    config["datadir"]=data_dir
