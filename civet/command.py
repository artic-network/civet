#!/usr/bin/env python3
from civet import __version__

import setuptools
import argparse
import os.path
import snakemake
import sys
import tempfile
import csv
import os
import yaml
from datetime import datetime
from Bio import SeqIO
import pkg_resources
from . import _program

from reportfunk.funks import io_functions as qcfunk
from reportfunk.funks import report_functions as rfunk
from reportfunk.funks import custom_logger as custom_logger
import civetfunks as cfunk

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='civet: Cluster Investivation & Virus Epidemiology Tool', 
    usage='''civet -i <input.csv> [options]
civet -i EDB3588,EDB3589 [options]
civet -i <config.yaml>
civet -fm adm2=Edinburgh sample_date=2020-03-01:2020-04-01 [options]''')

    io_group = parser.add_argument_group('input output options')
    io_group.add_argument('-i',"--input", action="store",help="Input config file in yaml format, csv file (with minimally an input_column header, Default=`name`) or comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="input")
    io_group.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the large phylogeny. E.g. -fm adm2=Edinburgh sample_date=2020-03-01:2020-04-01")
    io_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('-f','--fasta', action="store",help="Optional fasta query.", dest="fasta")
    io_group.add_argument('--max-ambig', action="store", type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="maxambig")
    io_group.add_argument('--min-length', action="store", type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="minlen")

    data_group = parser.add_argument_group('data source options')
    data_group.add_argument('-d','--datadir', action="store",help="Local directory that contains the data files. Default: civet-cat")
    data_group.add_argument('--CLIMB', action="store_true",dest="climb",help="Indicates you're running CIVET from within CLIMB, uses default paths in CLIMB to access data")
    data_group.add_argument("-r",'--remote-sync', action="store_true",dest="remote",help="Remotely access lineage trees from CLIMB")
    data_group.add_argument("-uun","--your-user-name", action="store", help="Your CLIMB COG-UK username. Required if running with --remote-sync flag", dest="uun")
    data_group.add_argument('--input-column', action="store",help="Column in input csv file to match with database. Default: name", dest="input_column")
    data_group.add_argument('--data-column', action="store",help="Option to search COG database for a different id type. Default: COG-UK ID", dest="data_column")
    data_group.add_argument("--display-name", action="store", help="Column in input csv file with display names for seqs. Default: same as input column", dest="display_name")

    report_group = parser.add_argument_group('report customisation')
    report_group.add_argument('-sc',"--sequencing-centre", action="store",help="Customise report with logos from sequencing centre.", dest="sequencing_centre")
    report_group.add_argument('--colour-by', action="store", help="Comma separated string of fields to display as coloured dots rather than text in report trees. Optionally add colour scheme eg adm1=viridis", dest="colour_by")
    report_group.add_argument('--tree-fields', action="store",help="Comma separated string of fields to display in the trees in the report. Default: country", dest="tree_fields")
    report_group.add_argument('--label-fields', action="store", help="Comma separated string of fields to add to tree report labels.", dest="label_fields")
    report_group.add_argument("--date-fields", action="store", help="Comma separated string of metadata headers containing date information.", dest="date_fields")
    report_group.add_argument("--node-summary", action="store", help="Column to summarise collapsed nodes by. Default = Global lineage", dest="node_summary")
    report_group.add_argument("--table-fields", action="store", help="Fields to include in the table produced in the report", dest="table_fields")
    report_group.add_argument("--include-snp-table", action="store", help="Include information about closest sequence in database in table. Default is False", dest="include_snp_table")
    report_group.add_argument('--no-snipit', action="store_true",help="Don't run snipit graph", dest="no_snipit")
    report_group.add_argument('--include-bars', action="store_true",help="Render barcharts in the output report", dest="include_bars")
    report_group.add_argument('--cog-report', action="store_true",help="Run summary cog report. Default: outbreak investigation",dest="cog_report")
    report_group.add_argument('--omit-appendix', action="store_true", help="Omit the appendix section. Default=False", dest="omit_appendix")
    report_group.add_argument('--private', action="store_true", help="remove adm2 references from background sequences. Default=False")

    tree_group = parser.add_argument_group('tree context options')
    tree_group.add_argument('--distance', action="store",help="Extraction from large tree radius. Default: 2", dest="distance",type=int)
    tree_group.add_argument('--up-distance', action="store",help="Upstream distance to extract from large tree. Default: 2", dest="up_distance",type=int)
    tree_group.add_argument('--down-distance', action="store",help="Downstream distance to extract from large tree. Default: 2", dest="down_distance",type=int)
    tree_group.add_argument('--collapse-threshold', action='store',type=int,help="Minimum number of nodes to collapse on. Default: 1", dest="threshold")

    map_group = parser.add_argument_group('map rendering options')
    map_group.add_argument('--local-lineages',action="store_true",dest="local_lineages",help="Contextualise the cluster lineages at local regional scale. Requires at least one adm2 value in query csv.")
    map_group.add_argument('--date-restriction',action="store_true",dest="date_restriction",help="Chose whether to date-restrict comparative sequences at regional-scale.")
    map_group.add_argument('--date-range-start',action="store", type=str, dest="date_range_start", help="Define the start date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    map_group.add_argument('--date-range-end', action="store", type=str, dest="date_range_end", help="Define the end date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    map_group.add_argument('--date-window',action="store", type=int, dest="date_window",help="Define the window +- either side of cluster sample collection date-range. Default is 7 days.")
    map_group.add_argument("--map-sequences", action="store_true", dest="map_sequences", help="Map the sequences themselves by adm2, coordinates or otuer postcode.")
    map_group.add_argument("--map-info", required=False, dest="map_info", help="columns containing EITHER x and y coordinates as a comma separated string OR outer postcodes for mapping sequences OR Adm2")
    map_group.add_argument("--input-crs", required=False, dest="input_crs", help="Coordinate reference system for sequence coordinates")
    map_group.add_argument("--colour-map-by", required=False, dest="colour_map_by", help="Column to colour mapped sequences by")
    
    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument('-b','--launch-browser', action="store_true",help="Optionally launch md viewer in the browser using grip",dest="launch_browser")
    misc_group.add_argument('-c','--generate-config',dest="generate_config",action="store_true",help="Rather than running a civet report, generate a config file based on the command line arguments provided")
    misc_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    misc_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads")
    misc_group.add_argument("-v","--version", action='version', version=f"civet {__version__}")
    
    """
    Exit with help menu if no args supplied
    """
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    """
    Initialising dicts
    """
    # create the config dict to pass through to the snakemake file
    config = {}
    # get the default values from civetfunks
    default_dict = cfunk.get_defaults()

    """
    Output directory 
    (needed prior to -i, in case input ids need to be written to a file)
    """
    # default output dir
    qcfunk.get_outdir(args.outdir,cwd,config)

    """
    Input file (-i/--input) 
    Valid inputs are input.csv; ID1,ID2,ID3; config.yaml/config.yml
    
    If there's an input fasta file- add to the config dict

    """
    # find the query csv, or string of ids, or config file
    query,configfile = qcfunk.type_input_file(args.input,cwd,config)

    # if a yaml file is detected, add everything in it to the config dict
    if configfile:
        config = qcfunk.parse_yaml_file(configfile, config)
    
    
    """
    Get tempdir and data dir. 
    Check if data has the right columns needed.
    The following rely on things that come out of the 
    input config or csv files so shouldnt be moved up above that.

    - tempdir
    - datadir
    """
    
    # specifying temp directory, outdir if no_temp (tempdir becomes working dir)
    tempdir = qcfunk.get_temp_dir(args.tempdir, args.no_temp,cwd,config)

    # find the data dir
    cfunk.get_datadir(args.climb,args.uun,args.datadir,args.remote,cwd,config,default_dict)

    # add data and input columns to config
    qcfunk.data_columns_to_config(args,config,default_dict)

    # check if metadata has the right columns, background_metadata_header added to config
    qcfunk.check_metadata_for_seach_columns(config,default_dict)

    """
    from metadata parsing 

    relies on the data dir being found 
    """
    # generate query from metadata
    if args.from_metadata or "from_metadata" in config:
        metadata = config["metadata"]
        query = qcfunk.generate_query_from_metadata(args.from_metadata,metadata,config)

    """
    The query file could have been from one of
    - input.csv
    - id string input, created csv
    - from_metadata generated query csv

    (all either specified in config or via command line)
    """
    # check query exists or add ids to temp query file
    qcfunk.check_query_file(query, cwd, config)

    """
    Input fasta file 
    sourcing and qc checks
    """
    # find the query fasta
    qcfunk.get_query_fasta(args.fasta,cwd, config)
    
    # run qc on the input sequence file
    qcfunk.input_file_qc(args.minlen,args.maxambig,config,default_dict)

    """
    Accessing the civet package data and 
    selecting the mapping files, 
    the sequencing centre header
    """
    # accessing package data and adding to config dict
    cfunk.get_package_data(args.cog_report,thisdir,config,default_dict)

    """
    Report options and args added to config, seq header file retrieved
    """
    # check args, config, defaultdict for report group options
    cfunk.report_group_to_config(args,config,default_dict)

    # get seq centre header file from pkg data
    cfunk.get_sequencing_centre_header(config)

    
    """
    Mapping options parsing and 
    qc of the input
    """
    
    # check args, config, defaultdict for mapping group options
    cfunk.map_group_to_config(args,config,default_dict)

    # check args, config, defaultdict for data group options
    qcfunk.data_columns_to_config(args,config,default_dict)

    # parse the input csv, check col headers and get fields if fields specified
    qcfunk.check_label_and_tree_and_date_fields(args.tree_fields, args.label_fields,args.colour_by, args.date_fields, args.input_column, args.display_name, config,default_dict)
        
    # map sequences configuration
    cfunk.map_sequences_config(args.map_sequences,args.colour_map_by,args.map_info,args.input_crs,config,default_dict)
    
    # local lineages qc
    cfunk.local_lineages_qc(config,default_dict)

    """
    Parsing the tree_group arguments, 
    config or default options
    """

    # global now the only search option
    cfunk.define_seq_db(config,default_dict)

    # extraction radius configuration
    qcfunk.distance_config(args.distance,args.up_distance,args.down_distance,config,default_dict) 

    """
    Parsing the report_group arguments, 
    config or default options
    """

    # make title
    rfunk.make_title(config)
    # deal with free text
    rfunk.free_text_args(config)
    rfunk.appendix(args.omit_appendix,config)

    #deal with figures
    rfunk.bars(args.include_bars, config)

    #get table headers
    qcfunk.check_table_fields(args.table_fields, args.include_snp_table, config,default_dict)
        
    # summarising collapsed nodes config
    qcfunk.node_summary(args.node_summary,config)

    for key in default_dict:
        if key not in config:
            config[key] = default_dict[key]

    """
    Miscellaneous options parsing

    """
    if args.launch_browser:
        config["launch_browser"]=True

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]=" "
    else:
        quiet_mode = True
        config["quiet_mode"]="--quiet"

    threads = qcfunk.check_arg_config_default("threads",args.threads,config,default_dict)
    config["threads"]= int(threads)

    if args.generate_config:
        qcfunk.make_config_file("civet_config.yaml",config)

    
    # find the master Snakefile
    snakefile = qcfunk.get_snakefile(thisdir)

    if args.verbose:
        
        for k in sorted(config):
            print(qcfunk.green(k), config[k])
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=tempdir,config=config, cores=threads,lock=False
                                    )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
                                    config=config, cores=threads,lock=False,quiet=True,log_handler=logger.log_handler
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()