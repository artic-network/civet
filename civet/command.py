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
from datetime import date
from Bio import SeqIO
import pkg_resources
from . import _program

from reportfunk.funks import io_functions as qcfunk
from reportfunk.funks import report_functions as rfunk
from reportfunk.funks import custom_logger as custom_logger
from reportfunk.funks import log_handler_handle as lh
import civetfunks as cfunk

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()
today = date.today()

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'


def get_global_snipit_snakefile(run_dir):
    snakefile = os.path.join(run_dir, 'scripts','global_snipit.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n')
        sys.exit(-1)
    return snakefile


def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False, prog = _program, 
    description=cfunk.preamble(__version__), 
    usage='''
\tcivet -i <config.yaml> [options]
\tcivet -i input.csv [options]
\tcivet -i ID1,IS2 [options]
\tcivet -fm <column=match> [options]\n\n''')

    io_group = parser.add_argument_group('input output options')
    io_group.add_argument('-i',"--input", action="store",help="Input config file in yaml format, csv file (with minimally an input_column header, Default=`name`) or comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="input")
    io_group.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the large phylogeny. E.g. -fm adm2=Edinburgh sample_date=2020-03-01:2020-04-01")
    io_group.add_argument('-o','--output-prefix',action="store",help="Prefix of output directory & report name: Default: civet",dest="output_prefix")
    io_group.add_argument('--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('-f','--fasta', action="store",help="Optional fasta query.", dest="fasta")
    io_group.add_argument('--max-ambiguity', action="store", type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="max_ambiguity")
    io_group.add_argument('--min-length', action="store", type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="min_length")
    
    data_group = parser.add_argument_group('data source options')
    data_group.add_argument('-d','--datadir', action="store",help="Local directory that contains the data files. Default: civet-cat")
    data_group.add_argument("-m","--background-metadata",action="store",dest="background_metadata",help="Custom metadata file that corresponds to the large global tree/ alignment. Should have a column `sequence_name`.")
    data_group.add_argument("--background-tree", action="store", dest="background_tree", help="Custom tree file.")
    data_group.add_argument("--background-sequences", action="store", dest="background_sequences", help="Custom background fasta file.")
    data_group.add_argument('--CLIMB', action="store_true",dest="climb",help="Indicates you're running CIVET from within CLIMB, uses default paths in CLIMB to access data")
    data_group.add_argument("-r",'--remote', action="store_true",dest="remote",help="Remotely access lineage trees from CLIMB")
    data_group.add_argument("-uun","--your-user-name", action="store", help="Your CLIMB COG-UK username. Required if running with --remote flag", dest="uun")
    data_group.add_argument('--input-column', action="store",help="Column in input csv file to match with database. Default: name", dest="input_column")
    data_group.add_argument('--data-column', action="store",help="Option to search COG database for a different id type. Default: COG-UK ID", dest="data_column")

    report_group = parser.add_argument_group('report customisation')
    report_group.add_argument('-sc',"--sequencing-centre", action="store",help="Customise report with logos from sequencing centre.", dest="sequencing_centre")
    report_group.add_argument("--display-name", action="store", help="Column in input csv file with display names for seqs. Default: same as input column", dest="display_name")
    report_group.add_argument("--sample-date-column", action="store", help="Column in input csv with sampling date in it. Default='sample_date'", dest="sample_date_column")
    report_group.add_argument("--database-sample-date-column", action="store", help="Colum in background metadata containing sampling date. Default='sample_date'", dest="database_sample_date_column")
    report_group.add_argument('--colour-by', action="store", help="Comma separated string of fields to display as coloured dots rather than text in report trees. Optionally add colour scheme eg adm1=viridis", dest="colour_by")
    report_group.add_argument('--tree-fields', action="store",help="Comma separated string of fields to display in the trees in the report. Default: country", dest="tree_fields")
    report_group.add_argument('--label-fields', action="store", help="Comma separated string of fields to add to tree report labels.", dest="label_fields")
    report_group.add_argument("--date-fields", action="store", help="Comma separated string of metadata headers containing date information.", dest="date_fields")
    report_group.add_argument("--node-summary", action="store", help="Column to summarise collapsed nodes by. Default = Global lineage", dest="node_summary")
    report_group.add_argument("--table-fields", action="store", help="Fields to include in the table produced in the report. Query ID, name of sequence in tree and the local tree it's found in will always be shown", dest="table_fields")
    report_group.add_argument("--remove-snp-table", action="store_true", help="Include information about closest sequence in database in table. Default is False", dest="remove_snp_table")
    report_group.add_argument('--no-snipit', action="store_true",help="Don't run snipit graph", dest="no_snipit")
    report_group.add_argument('--global_snipit', action="store_true", help="Create a global snipit figure for all focal sequences", dest="global_snipit")
    report_group.add_argument('--include-bars', action="store_true",help="Render barcharts in the output report", dest="include_bars")
    report_group.add_argument('--omit-appendix', action="store_true", help="Omit the appendix section. Default=False", dest="omit_appendix")
    report_group.add_argument('--omit-trees', action="store_true", help="Omit trees.", dest="omit_trees")
    report_group.add_argument('--context-table-summary', help="Provide a field to summarise the context by", dest='context_table_summary')

    tree_group = parser.add_argument_group('tree context options')
    tree_group.add_argument('--distance', action="store",help="Extraction from large tree radius. Default: 2", dest="distance",type=int)
    tree_group.add_argument('--up-distance', action="store",help="Upstream distance to extract from large tree. Default: 2", dest="up_distance",type=int)
    tree_group.add_argument('--down-distance', action="store",help="Downstream distance to extract from large tree. Default: 2", dest="down_distance",type=int)
    tree_group.add_argument('--collapse-threshold', action='store',help="Minimum number of nodes to collapse on. Default: 1", dest="collapse_threshold",type=int)
    tree_group.add_argument('-p','--protect',nargs='*', dest="protect",help="Protect nodes from collapse if they match the search query in the metadata file supplied. E.g. -p adm2=Edinburgh sample_date=2020-03-01:2020-04-01")

    map_group = parser.add_argument_group('map rendering options')
    map_group.add_argument('--local-lineages',action="store_true",dest="local_lineages",help="Contextualise the cluster lineages at local regional scale. Requires at least one adm2 value in query csv.")
    map_group.add_argument('--date-restriction',action="store_true",dest="date_restriction",help="Chose whether to date-restrict comparative sequences at regional-scale.")
    map_group.add_argument('--date-range-start',action="store", type=str, dest="date_range_start", help="Define the start date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    map_group.add_argument('--date-range-end', action="store", type=str, dest="date_range_end", help="Define the end date from which sequences will COG sequences will be used for local context. YYYY-MM-DD format required.")
    map_group.add_argument('--date-window',action="store", type=int, dest="date_window",help="Define the window +- either side of cluster sample collection date-range. Default is 7 days.")
    map_group.add_argument("--map-sequences", action="store_true", dest="map_sequences", help="Map the sequences themselves by adm2, coordinates or outer postcode.")
    map_group.add_argument("--map-info", required=False, dest="map_info", help="columns containing EITHER x and y coordinates as a comma separated string OR outer postcodes for mapping sequences OR Adm2")
    map_group.add_argument("--input-crs", required=False, dest="input_crs", help="Coordinate reference system for sequence coordinates")
    map_group.add_argument("--colour-map-by", required=False, dest="colour_map_by", help="Column to colour mapped sequences by")
    
    run_group = parser.add_argument_group('run options')
    run_group.add_argument("--cluster",action="store_true",help="Run cluster civet pipeline. Requires -i/--input csv",dest="cluster")
    run_group.add_argument("--update",action="store_true",help="Check for changes from previous run of civet. Requires -fm/--from-metadata option in a config.yaml file from previous run",dest="update")
    run_group.add_argument("--udpate",action="store_true",help="Check for changes from previous run of civet. Requires -fm/--from-metadata option in a config.yaml file from previous run",dest="udpate")
    run_group.add_argument('-c','--generate-config',dest="generate_config",action="store_true",help="Rather than running a civet report, just generate a config file based on the command line arguments provided")
    run_group.add_argument('-b','--launch-browser', action="store_true",help="Optionally launch md viewer in the browser using grip",dest="launch_browser")

    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument("--safety-level", action="store", type=int, dest="safety_level",help="Level of anonymisation for users. Options: 0 (no anonymity), 1 (no COGIDs on background data), 2 (no adm2 on data). Default: 1")
    misc_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    misc_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.",dest="no_temp")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument("--art",action="store_true",help="Print art")
    misc_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads")
    misc_group.add_argument("-v","--version", action='version', version=f"civet {__version__}")
    misc_group.add_argument("-h","--help",action="store_true",dest="help")
    
    """
    Exit with help menu if no args supplied
    """
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)
    
    if args.art:
        cfunk.be_arty()
        sys.exit(0) 
    
    """
    Initialising dicts
    """

    # get the default values from civetfunk
    config = cfunk.get_defaults()

    """
    Input file (-i/--input) 
    Valid inputs are input.csv; ID1,ID2,ID3; config.yaml/config.yml
    
    If there's an input fasta file- add to the config dict

    """
    # find the query csv, or string of ids, or config file
    query,configfile = qcfunk.type_input_file(args.input,cwd,config)

    # if a yaml file is detected, add everything in it to the config dict
    if configfile:
        qcfunk.parse_yaml_file(configfile, config)
    
    """
    Report options and args added to config, seq header file retrieved
    """
    # check args for report group options
    cfunk.report_group_to_config(args,config)


    # update and cluster options

    cfunk.configure_update(args.update,args.udpate,config)
    qcfunk.add_arg_to_config("cluster",args.cluster, config)

    """
    Get outdir, tempdir and data dir. 
    Check if data has the right columns needed.
    The following rely on things that come out of the 
    input config or csv files so shouldnt be moved up above that.

    - tempdir
    - datadir
    """
    # default output dir
    qcfunk.get_outdir(args.outdir,args.output_prefix,cwd,config)

    # specifying temp directory, outdir if no_temp (tempdir becomes working dir)
    tempdir = qcfunk.get_temp_dir(args.tempdir, args.no_temp,cwd,config)

    qcfunk.add_arg_to_config("remote",args.remote, config)

    # find the data dir
    cfunk.get_datadir(args.climb,args.uun,args.datadir,args.background_metadata, args.background_tree, args.background_sequences, cwd,config)

    # add data and input columns to config
    qcfunk.data_columns_to_config(args,config)

    # check if metadata has the right columns, background_metadata_header added to config
    qcfunk.check_metadata_for_search_columns(config)


    """
    from metadata parsing 

    relies on the data dir being found 
    """
    # generate query from metadata
    qcfunk.add_arg_to_config("from_metadata",args.from_metadata, config)
    if config["from_metadata"]:
        
        qcfunk.from_metadata_checks(config)

        metadata = config["background_metadata"]
        config["no_snipit"]=True

        if config["update"]:
            query_file = os.path.join(config["outdir"], "update_query.csv")
            run_update = cfunk.check_for_update(query_file,config)
            if not run_update:
                print(qcfunk.cyan('Note: no new sequences to report.\nExiting.'))
                sys.exit(0)
            else:
                query = config["query"] # gets added updated in the check_for_update function
        else:
            query_file = os.path.join(config["outdir"], "from_metadata_query.csv")
            query = qcfunk.generate_query_from_metadata(query_file,args.from_metadata,metadata,config)
    else:
        if config["update"]:
            cfunk.check_update_dependencies(config)
            
    """
    The query file could have been from one of
    - input.csv
    - id string input, created csv
    - from_metadata generated query csv

    (all either specified in config or via command line)
    """
    # check query exists or add ids to temp query file
    qcfunk.check_query_file(query, cwd, config)

    # check if metadata has the right columns, background_metadata_header added to config
    qcfunk.check_query_for_input_column(config)

    """
    Input fasta file 
    sourcing and qc checks
    """
    # find the query fasta
    qcfunk.get_query_fasta(args.fasta,cwd, config)
    
    # run qc on the input sequence file
    num_seqs = qcfunk.input_file_qc(args.min_length,args.max_ambiguity,config)
    
    """
    Quick check in background data
    """
    if num_seqs == 0:
        # check if any queries in background or if fasta supplied
        qcfunk.check_background_for_queries(config)

    """
    Accessing the civet package data and 
    selecting the mapping files, 
    the sequencing centre header
    """
    # accessing package data and adding to config dict
    cfunk.get_package_data(thisdir,config)

    # get seq centre header file from pkg data
    cfunk.get_sequencing_centre_header(config)

    
    """
    Mapping options parsing and 
    qc of the input
    """
    
    # check args for mapping group options
    cfunk.map_group_to_config(args,config)

    # check args for data group options
    qcfunk.data_columns_to_config(args,config)

    # parse the input csv, check col headers and get fields if fields specified
    qcfunk.check_label_and_tree_and_date_fields(config)
        
    # map sequences configuration
    cfunk.map_sequences_config(config)
    
    # local lineages qc
    cfunk.local_lineages_qc(config)

    #check adm2s
    if config["local_lineages"] or config["map_sequences"]:
        cfunk.check_adm2_values(config)

    """
    Parsing the tree_group arguments, 
    config or default options
    """

    # global now the only search option
    cfunk.define_seq_db(config)

    # extraction radius configuration
    qcfunk.distance_config(args.distance,args.up_distance,args.down_distance,config) 

    # extraction radius configuration
    qcfunk.collapse_config(args.collapse_threshold,config) 


    qcfunk.parse_protect(args.protect,config["background_metadata"],config)

    """
    Parsing the report_group arguments, 
    config or default options
    """

    # make title
    rfunk.make_title(config)
    # deal with free text
    rfunk.free_text_args(config)

    #get table headers
    qcfunk.check_table_fields(config) 
        
    # summarising collapsed nodes config
    qcfunk.check_summary_field("node_summary",config)

    qcfunk.collapse_summary_path_to_config(config)

   

    """
    Miscellaneous options parsing

    """

    qcfunk.add_arg_to_config("launch_browser",args.launch_browser,config)

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["log_string"] = ""
    else:
        quiet_mode = True
        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "

    qcfunk.add_arg_to_config("threads",args.threads,config)
    
    try:
        config["threads"]= int(config["threads"])
    except:
        sys.stderr.write(qcfunk.cyan('Error: Please specifiy an integer for variable `threads`.\n'))
        sys.exit(-1)
    threads = config["threads"]


    if args.safety_level != None:
        config["safety_level"]= args.safety_level
    
    if config["remote"]:
        config["safety_level"] = 2
    
    try:
        safety_level = int(config["safety_level"])
    except:
        sys.stderr.write(qcfunk.cyan('Error: Please specifiy either 0, 1 or 2 for variable `safety_level`.\n'))
        sys.exit(-1)

    if safety_level in [0,1,2]:
        config["safety_level"]= int(config["safety_level"])
    else:
        sys.stderr.write(qcfunk.cyan('Error: Please specifiy either 0, 1 or 2 for variable `safety_level`.\n'))
        sys.exit(-1)

    if args.generate_config:
        qcfunk.make_config_file("civet_config.yaml",config)

    if config["global_snipit"]:
        qcfunk.add_arg_to_config("global_snipit", args.global_snipit, config)
        global_snipit_snakefile = get_global_snipit_snakefile(thisdir)
        status = snakemake.snakemake(global_snipit_snakefile, printshellcmds=True, forceall=True,
                                     force_incomplete=True,

                                     workdir=tempdir, config=config, cores=threads, lock=False)

        if not status:
            print(qcfunk.cyan(f"Error: Cluster civet did not successfully run"))
            sys.exit(-1)

    """
    cluster civet checks
    - arg, config, default
    - is there a query?
    - check if new things in the local tree
    - if new sequences, run main civet with updated query
    - if no new sequences, exit
    
    """
    # cluster civet 

    if config["cluster"]:

        config["today"] = today
        cfunk.configure_cluster(config)

        cluster_snakefile = qcfunk.get_cluster_snakefile(thisdir)


        if args.verbose:
            print("\n**** CONFIG ****")
            for k in sorted(config):
                print(qcfunk.green(k), config[k])
            status = snakemake.snakemake(cluster_snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        workdir=tempdir,config=config, cores=threads,lock=False
                                        )
        else:
            logger = custom_logger.Logger()
            status = snakemake.snakemake(cluster_snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
                                        config=config, cores=threads,lock=False,quiet=True,log_handler=logger.log_handler
                                        )

        if not status:
            print(qcfunk.cyan(f"Error: Cluster civet did not successfully run"))
            sys.exit(-1)

        new_seqs, cluster_csv = cfunk.check_for_new_in_cluster(config)

        print(qcfunk.green(f"\nNew sequences found in cluster {today}: ") + f"{new_seqs}")

        if not new_seqs:
            print(qcfunk.cyan(f"No new sequences in cluster today, {today}"))
            sys.exit(0)
        else:
            config["query"] = cluster_csv

    # find the master Snakefile
    snakefile = qcfunk.get_snakefile(thisdir)

    if args.verbose:
        print("\n**** CONFIG ****")
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