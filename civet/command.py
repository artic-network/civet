#!/usr/bin/env python3
from civet import __version__
import setuptools
import argparse
import os.path
import snakemake
import sys
from tempfile import gettempdir
import tempfile
import pprint
import json
import csv
import input_qc_functions as qcfunk
import os
from datetime import datetime
from Bio import SeqIO

import pkg_resources
from . import _program

"""
Need query_csv, metadata, fasta (opt)
"""


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='civet: Cluster Investivation & Virus Epidemiology Tool', 
    usage='''civet <query> [options]''')

    parser.add_argument('query',help="Input csv file with minimally `name` as a column header. Can include additional fields to be incorporated into the analysis, e.g. `sample_date`")
    parser.add_argument('-i',"--id-string", action="store_true",help="Indicates the input is a comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="ids")
    parser.add_argument('-f','--fasta', action="store",help="Optional fasta query.", dest="fasta")
    parser.add_argument('-sc',"--sequencing-centre", action="store",help="Customise report with logos from sequencing centre.", dest="sequencing_centre")
    parser.add_argument('--CLIMB', action="store_true",dest="climb",help="Indicates you're running CIVET from within CLIMB, uses default paths in CLIMB to access data")
    parser.add_argument('--cog-report', action="store_true",help="Run summary cog report. Default: outbreak investigation",dest="cog_report")
    parser.add_argument("-r",'--remote-sync', action="store_true",dest="remote",help="Remotely access lineage trees from CLIMB")
    parser.add_argument("-uun","--your-user-name", action="store", help="Your CLIMB COG-UK username. Required if running with --remote-sync flag", dest="uun")
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('-b','--launch-browser', action="store_true",help="Optionally launch md viewer in the browser using grip",dest="launch_browser")
    parser.add_argument('-d','--datadir', action="store",help="Local directory that contains the data files")
    parser.add_argument('--fields', action="store",help="Comma separated string of fields to display in the trees in the report. Default: country")
    parser.add_argument('--display', action="store", help="Comma separated string of fields to display as coloured dots rather than text in report trees. Optionally add colour scheme eg adm1=viridis", dest="display")
    parser.add_argument('--label-fields', action="store", help="Comma separated string of fields to add to tree report labels.", dest="label_fields")
    parser.add_argument("--node-summary", action="store", help="Column to summarise collapsed nodes by. Default = Global lineage", dest="node_summary")
    parser.add_argument('--input-column', action="store",help="Column in input csv file to match with database. Default: name", dest="input_column",default="name")
    parser.add_argument('--search-field', action="store",help="Option to search COG database for a different id type. Default: COG-UK ID", dest="data_column",default="central_sample_id")
    parser.add_argument('--distance', action="store",help="Extraction from large tree radius. Default: 2", dest="distance",default=2)
    parser.add_argument('--up-distance', action="store",help="Upstream distance to extract from large tree. Default: 2", dest="up_distance",default=2)
    parser.add_argument('--down-distance', action="store",help="Downstream distance to extract from large tree. Default: 2", dest="down_distance",default=2)
    parser.add_argument('--add-bars', action="store_true",help="Render boxplots in the output report", dest="add_bars")
    # parser.add_argument('--delay-tree-collapse',action="store_true",dest="delay_tree_collapse",help="Wait until after iqtree runs to collapse the polytomies. NOTE: This may result in large trees that take quite a while to run.")
    parser.add_argument('-g','--global',action="store_true",dest="search_global",help="Search globally.")
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    parser.add_argument('--collapse-threshold', action='store',type=int,help="Minimum number of nodes to collapse on. Default: 1", dest="threshold", default=1)
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument('--max-ambig', action="store", default=0.5, type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=10000, type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="minlen")
    parser.add_argument('--local-lineages',action="store_true",dest="local_lineages",help="Contextualise the cluster lineages at local regional scale. Requires at least one adm2 value in query csv.", default=False)
    parser.add_argument('--date-restriction',action="store_true",dest="date_restriction",help="Chose whether to date-restrict comparative sequences at regional-scale.", default=False)
    parser.add_argument('--date-range-start',action="store",default="None", type=str, dest="date_range_start", help="Define the start date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    parser.add_argument('--date-range-end', action="store", default="None", type=str, dest="date_range_end", help="Define the end date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    parser.add_argument('--date-window',action="store",default=7, type=int, dest="date_window",help="Define the window +- either side of cluster sample collection date-range. Default is 7 days.")
    parser.add_argument("-v","--version", action='version', version=f"civet {__version__}")

    parser.add_argument("--map-sequences", action="store_true", dest="map_sequences", help="Map the coordinate points of sequences, coloured by a triat.")
    parser.add_argument("--x-col", required=False, dest="x_col", help="column containing x coordinate for mapping sequences")
    parser.add_argument("--y-col", required=False, dest="y_col", help="column containing y coordinate for mapping sequences")
    parser.add_argument("--input-crs", required=False, dest="input_crs", help="Coordinate reference system of sequence coordinates")
    parser.add_argument("--mapping-trait", required=False, dest="mapping_trait", help="Column to colour mapped sequences by")
    
    # Exit with help menu if no args supplied
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    # create the config dict to pass through to the snakemake file
    config = {
        "trim_start":265,   # where to pad to using datafunk
        "trim_end":29674,   # where to pad after using datafunk
        "search_field":args.data_column,
        "input_column":args.input_column,
        "force":"True",
        "date_range_start":args.date_range_start,
        "date_range_end":args.date_range_end,
        "date_window":args.date_window
        }
    # find the query csv, or string of ids, or config file
    query = qcfunk.parse_input_query(args.query,args.ids,cwd,config)
    
    # find the master Snakefile
    snakefile = qcfunk.get_snakefile(thisdir)
    
    # find the query fasta
    qcfunk.get_query_fasta(args.fasta,cwd, config)

    # default output dir
    outdir = qcfunk.get_outdir(args.outdir,cwd,config)

    # specifying temp directory, outdir if no_temp
    tempdir =qcfunk.get_temp_dir(args.tempdir, args.no_temp,cwd,config)

    # parse the input csv, check col headers and get fields if fields specified

    qcfunk.check_label_and_colour_fields(query, args.query, args.fields, args.label_fields,args.display, args.input_column, config)
        
    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print(f"Number of threads: {threads}\n")
    
    if args.add_bars:
        config["add_bars"]= True
    else:
        config["add_bars"]= False
        
    qcfunk.map_sequences_config(args.map_sequences,args.mapping_trait,args.x_col,args.y_col,args.input_crs,query,config)
    
    qcfunk.local_lineages_config(args.local_lineages,query,config)
    
    if args.date_restriction:
        config['date_restriction'] = "True"
    else:
        config['date_restriction'] = "False"

    """ 
    QC steps:
    1) check csv header
    2) check fasta file N content
    3) write a file that contains just the seqs to run
    """

    # find the data dir
    data_dir = qcfunk.get_datadir(args.climb,args.datadir,args.remote,cwd,config)
    # if remote flag, and uun provided, sync data from climb
    if args.remote:
        config["remote"]= "True"
        if args.uun:
            config["username"] = args.uun
        
            rsync_command = f"rsync -avzh {args.uun}@bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{data_dir}'"
            print(f"Syncing civet data to {data_dir}")
            status = os.system(rsync_command)
            if status != 0:
                sys.stderr.write("Error: rsync command failed.\nCheck your user name is a valid CLIMB username e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and are in the UK\n\n")
                sys.exit(-1)
        else:
            rsync_command = f"rsync -avzh bham.covid19.climb.ac.uk:/cephfs/covid/bham/civet-cat '{data_dir}'"
            print(f"Syncing civet data to {data_dir}")
            status = os.system(rsync_command)
            if status != 0:
                sys.stderr.write("Error: rsync command failed.\nCheck your ssh is configured with Host bham.covid19.climb.ac.uk\nAlternatively enter your CLIMB username with -uun e.g. climb-covid19-smithj\nAlso, check if you have access to CLIMB from this machine and are in the UK\n\n")
                sys.exit(-1)
        cog_metadata,all_cog_metadata,cog_global_metadata = ("","","")
        cog_seqs,all_cog_seqs = ("","")
        cog_tree = ""

        cog_seqs = os.path.join(data_dir,"civet-cat","cog_alignment.fasta")
        all_cog_seqs = os.path.join(data_dir,"civet-cat","cog_alignment_all.fasta")

        cog_metadata = os.path.join(data_dir,"civet-cat","cog_metadata.csv")
        all_cog_metadata = os.path.join(data_dir,"civet-cat","cog_metadata_all.csv")

        cog_global_metadata = os.path.join(data_dir,"civet-cat","cog_global_metadata.csv")
        cog_global_seqs= os.path.join(data_dir,"civet-cat","cog_global_alignment.fasta")

        cog_tree = os.path.join(data_dir,"civet-cat","cog_global_tree.nexus")

        config["cog_seqs"] = cog_seqs
        config["all_cog_seqs"] = all_cog_seqs

        config["cog_metadata"] = cog_metadata
        config["all_cog_metadata"] = all_cog_metadata
        config["cog_global_metadata"] = cog_global_metadata
        config["cog_global_seqs"] = cog_global_seqs
        config["cog_tree"] = cog_tree

        print("Found cog data:")
        print("    -",cog_seqs)
        print("    -",all_cog_seqs)
        print("    -",cog_metadata)
        print("    -",all_cog_metadata)
        print("    -",cog_global_metadata)
        print("    -",cog_tree,"\n")

    elif not args.datadir and not args.climb:
        sys.stderr.write("""Error: no way to find source data.\n\nTo run civet please either\n1) ssh into CLIMB and run with --CLIMB flag\n\
2) Run using `--remote-sync` flag and your CLIMB username specified e.g. `-uun climb-covid19-otoolexyz`\n\
3) Specify a local directory with the appropriate files on. The following files are required:\n\
- cog_global_tree.nexus\n\
- cog_metadata.csv\n\
- cog_metadata_all.csv\n\
- cog_global_metadata.csv\n\
- cog_global_alignment.fasta\n\
- cog_alignment.fasta\n\n""")
        sys.exit(-1)


    # run qc on the input sequence file
    qcfunk.input_file_qc(args.fasta,args.minlen,args.maxambig,config)
    
    if args.search_global:
        config["global"]="True"
    else:
        config["global"]="False"

    config["delay_collapse"] = False

    # accessing package data and adding to config dict
    qcfunk.get_package_data(args.cog_report,thisdir,config)

    with open(cog_global_metadata, newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames

        if not args.node_summary:
                summary = "country"
        else:
            if args.node_summary in column_names:
                summary = args.node_summary
            else:
                sys.stderr.write(f"Error: {args.node_summary} field not found in metadata file\n")
                sys.exit(-1)
        
        print(f"Going to summarise collapsed nodes by: {summary}")
        config["node_summary"] = summary

    

    sc_list = ["PHEC", 'LIVE', 'BIRM', 'PHWC', 'CAMB', 'NORW', 'GLAS', 'EDIN', 'SHEF',
                 'EXET', 'NOTT', 'PORT', 'OXON', 'NORT', 'NIRE', 'GSTT', 'LOND', 'SANG',"NIRE"]

    if args.sequencing_centre:
        if args.sequencing_centre in sc_list:
            relative_file = os.path.join("data","headers",f"{args.sequencing_centre}.png")
            header = pkg_resources.resource_filename('civet', relative_file)
            print(f"using header file from {header}\n")
            config["sequencing_centre"] = header
        else:
            sc_string = "\n".join(sc_list)
            sys.stderr.write(f'Error: sequencing centre must be one of the following:\n{sc_string}\n')
            sys.exit(-1)
    else:
        relative_file = os.path.join("data","headers","DEFAULT.png")
        header = pkg_resources.resource_filename('civet', relative_file)
        print(f"using header file from {header}\n")
        config["sequencing_centre"] = header

    if args.distance:
        try:
            distance = int(args.distance) 
            config["up_distance"] = args.distance
            config["down_distance"] = args.distance
        except:
            sys.stderr.write('Error: distance must be an integer\n')
            sys.exit(-1)
    else:
        config["up_distance"] = "1"
        config["down_distance"] = "1"

    if args.up_distance:
        try:
            distance = int(args.up_distance) 
            config["up_distance"] = args.up_distance
        except:
            sys.stderr.write('Error: up-distance must be an integer\n')
            sys.exit(-1)

    if args.down_distance:
        try:
            distance = int(args.down_distance) 
            config["down_distance"] = args.down_distance
        except:
            sys.stderr.write('Error: down-distance must be an integer\n')
            sys.exit(-1)


    if args.threshold:
        try:
            threshold = int(args.threshold)
            config["threshold"] = args.threshold
        except:
            sys.stderr.write('Error: threshold must be an integer\n')
            sys.exit(-1)
    else:
        config["threshold"] = "1"

    if args.launch_browser:
        config["launch_browser"]="True"

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]="False"
    else:
        quiet_mode = True
        config["quiet_mode"]="True"

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=True,force_incomplete=True,workdir=tempdir,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()