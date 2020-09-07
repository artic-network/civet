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
    parser.add_argument('--distance', action="store",help="Extraction from large tree radius. Default: 2", dest="distance",type=int,default=2)
    parser.add_argument('--up-distance', action="store",help="Upstream distance to extract from large tree. Default: 2", dest="up_distance",type=int,default=2)
    parser.add_argument('--down-distance', action="store",help="Downstream distance to extract from large tree. Default: 2", dest="down_distance",type=int,default=2)
    parser.add_argument('--add-bars', action="store_true",help="Render boxplots in the output report", dest="add_bars")
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
        "date_window":args.date_window,
        "threshold": args.threshold
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
        
    # map sequences configuration
    qcfunk.map_sequences_config(args.map_sequences,args.mapping_trait,args.x_col,args.y_col,args.input_crs,query,config)
    
    # local lineages configuration
    qcfunk.local_lineages_config(args.local_lineages,query,config)
    
    if args.date_restriction:
        config['date_restriction'] = True
    else:
        config['date_restriction'] = False


    # find the data dir
    data_dir = qcfunk.get_datadir(args.climb,args.datadir,args.remote,cwd,config)
    # if remote flag, and uun provided, sync data from climb
    qcfunk.get_remote_data(args.remote,args.uun,data_dir,args.datadir,args.climb,config)

    # run qc on the input sequence file
    qcfunk.input_file_qc(args.fasta,args.minlen,args.maxambig,config)
    
    if args.search_global:
        config["global"]=True
    else:
        config["global"]=False

    config["delay_collapse"] = False

    # accessing package data and adding to config dict
    qcfunk.get_package_data(args.cog_report,thisdir,config)

    # summarising collapsed nodes config
    qcfunk.node_summary(args.node_summary,config)

    # get seq centre header file from pkg data
    qcfunk.get_sequencing_centre_header(args.sequencing_centre,config)
    
    # extraction radius configuration
    qcfunk.distance_config(args.distance, args.up_distance, args.down_distance, config)
 
    if args.launch_browser:
        config["launch_browser"]=True

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]=False
    else:
        quiet_mode = True
        config["quiet_mode"]=True

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=True,force_incomplete=True,workdir=tempdir,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()