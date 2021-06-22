#!/usr/bin/env python3
from civet import __version__

from civet.input_parsing import initialising as init
from civet.input_parsing import input_arg_parsing
from civet.input_parsing import data_arg_parsing
from civet.input_parsing import analysis_arg_parsing
from civet.input_parsing import input_data_parsing
from civet.input_parsing import report_arg_parsing

from civet.output_options import directory_setup
from civet.report_functions import global_report_functions

from civet.utils import misc
from civet.utils import dependency_checks
from civet.utils import data_install_checks

import civet.utils.custom_logger as custom_logger
from civet.utils.log_colours import green,cyan,red

import os
import sys
import argparse
import snakemake

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(add_help=False,
    description=misc.preamble(__version__), 
    usage='''
\tcivet -c <config.yaml> [options]
\tcivet -i input.csv [options]
\tcivet -ids ID1,ID2 [options]
\tcivet -fm <column=match> [options]\n\n''')

    i_group = parser.add_argument_group('Input options')
    i_group.add_argument('-c',"--config", action="store",help="Input config file in yaml format, all command line arguments can be passed via the config file.", dest="config")
    i_group.add_argument('-i',"--input-csv", action="store",help="Input csv file (with minimally an input_column header, Default=`name`)", dest="input_csv")
    i_group.add_argument('-icol',"--input-column", action="store",help="Column in input csv file to match with database. Default: `name`", dest="input_column")

    i_group.add_argument('-ids',"--id-string", action="store",help="Comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`.", dest="ids")

    i_group.add_argument('-f','--fasta', action="store",help="Optional fasta file. Sequence IDs must match to a query ID specified either in the input csv or ID string.", dest="fasta")
    i_group.add_argument('-n','--max-ambiguity', action="store", type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: 0.5",dest="max_ambiguity")
    i_group.add_argument('-l','--min-length', action="store", type=int,help="Minimum query length allowed to attempt analysis. Default: 20000",dest="min_length")
    
    i_group.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the large phylogeny. E.g. -fm adm2=Edinburgh sample_date=2020-03-01:2020-04-01")

    d_group = parser.add_argument_group('Background data options')
    d_group.add_argument('-d','--datadir', action="store",help="Directory containing the background data files.")
    d_group.add_argument("-bc","--background-csv",action="store",dest="background_csv",help="Custom metadata file for all background data. Should have a column matching '--data-column', Default: sequence_name.")
    d_group.add_argument("-bSNP","--background-SNPs",action="store",dest="background_SNPs",help="Optional SNP file for all background data. Civet will calculate this file if not supplied, which may take some time. Should have a column matching '--data-column', Default: sequence_name.")
    d_group.add_argument("-bf","--background-fasta", action="store", dest="background_fasta", help="Custom background fasta file for all background data. Sequence IDs should match the background metadata background_column.")
    d_group.add_argument("-bt","--background-tree", action="store", dest="background_tree", help="Custom background tree file for all background data. Tip names should match the background metadata background_column.")
    d_group.add_argument("-bcol",'--background-column', action="store",help="Column in background data to match with input IDs. Default: sequence_name", dest="background_column")
    d_group.add_argument("-fcol",'--fasta-column', action="store",help="Column in background data to match with input IDs. Default: `-dcol/--data-column`.", dest="fasta_column")

    o_group = parser.add_argument_group('Output options')
    o_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: civet-2021-XX-YY")
    o_group.add_argument('-p','--output-prefix',action="store",help="Prefix of output directory & report name: Default: civet",dest="output_prefix")
    o_group.add_argument('-ds','--datestamp', action="store",help="Append datestamp to directory name when using `-o/--outdir`. Default: `-o/--outdir` without a datestamp.")
    o_group.add_argument('--overwrite', action="store_true",help="Overwrite output directory. Default: append a number if `-o/--outdir` exists.")
    o_group.add_argument('--output-data',action="store_true",help="Output intermediate data files to the output directory",dest="output_data")
    o_group.add_argument('-temp','--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    o_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files. For development/ debugging purposes.",dest="no_temp")

    a_group = parser.add_argument_group("Analysis group")
    a_group.add_argument('-ts','--trim-start', type=int, action="store",dest="trim_start",help="Genome position to trim and pad to when aligning input sequences. Default: 265")
    a_group.add_argument('-te','--trim-end', type=int, action="store",dest="trim_end",help="Genome position to trim and pad from when aligning input sequences. Default: 29674")
    a_group.add_argument("-r","--reference-fasta",action="store",dest="reference_fasta",help="Custom reference genome to map and pad against. Must match the reference the background data was generated from.")
    a_group.add_argument('-ql','--query-limit', type=int, action="store",dest="query_limit",help="Max number of queries. Default: 5000")
    a_group.add_argument('-cs','--catchment-size', type=int, action="store",dest="catchment_size",help="Max number of sequences in a catchment. Default: 1000")
    a_group.add_argument('--downsample', type=int, action="store",dest="downsample",help="Indicates how many shared number of SNPs to downsample catchment by. Default: 2. Set to 0 to turn downsampling off.")

    r_group = parser.add_argument_group("Report options")
    r_group.add_argument("-rc", "--report-content", nargs='*', action="store", dest="report_content", help="""One or more comma separated numeric strings to define the report content. Default: 1,2,3""")
    r_group.add_argument("--anonymise", action="store_true", dest="anonymise",help="Generates arbitrary labels for sequences for dissemination")
    r_group.add_argument("-rcol", "--report-column", action="store", dest="report_column", help="Column containing alternative sequence names, for example patient IDs")
    r_group.add_argument("--table-content", action='store', dest="table_content", help="Columns to include in the table for queries. Default:--background_column,--background_date_column,source,lineage,country,catchment")
    r_group.add_argument("--catchment-table", action='store', dest="catchment_table", help="Columns to include in the summary table for catchments. Default: count,country,lineage")
    r_group.add_argument("-td", "--timeline-dates", action='store', dest="timeline_dates", help="Data to generate a timeline as a comma separated string")
    r_group.add_argument("--timeline-colours", action='store', dest="timeline_colours",help="Comma separated string of hex codes or names of HTML compatible colours to colour dots by in the timeline.")
    r_group.add_argument("-dcol","--date-column", action="store", dest="date_column", help="Column in input query with date data in. Default=sample_date")
    r_group.add_argument("-bdate", "--background-date-column", action="store", dest="background_date_column", help="Column in input data with date data in. Default=sample_date")
    r_group.add_argument("-lcol","--location-column", dest="location_column", action="store", help="Column in the background csv that contains geographical data. Default=country")


    m_group = parser.add_argument_group("Map options")
    m_group.add_argument("--map-file", action="store", dest="map_file", help="JSON or GeoJSON containing polygons to plot queries or background on.")
    m_group.add_argument("-lat","--latitude-column", dest="latitude_column", action="store", help="Column containing latitude coordinate information to plot queries on a map")
    m_group.add_argument("-long","--longitude-column", dest="longitude_column", action="store", help="Column containing longitude coordinate information to plot queries on a map")
    m_group.add_argument("-daterestric","--background-map-date-restriction", dest="background_map_date_restriction", action="store", help="Restrict the time frame for mapping background lineage diversity. Can be an integer (number of days either side of queries to restrict to) or a date range, format='YYYY-MM-DD:YYYY-MM-DD'")
    # m_group.add_argument("-maploc","--map-location", dest="map_location", action="store", help="Column in input metadata containing location to match to background  to map background lineage diversity by")
    m_group.add_argument("-maploc","--background-map-location", dest="background_map_location", action="store", help="Column in background metadata containing location to map background lineage diversity by")


    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument("--civet-mode", action="store", dest='civet_mode', help="if CLIMB then import UK specific modules. Default=Global")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument("--debug",action="store_true",help="Debugging mode.")
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

    dependency_checks.check_dependencies()
    
    # Initialise config dict
    config = init.setup_config_dict(cwd,args.config)
    
    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    snakefile = data_install_checks.get_snakefile(thisdir)
    
    # Threads and verbosity to config
    init.misc_args_to_config(args.verbose,args.threads,args.civet_mode,config)
    init.set_up_verbosity(config)

    # Analysis options, including ref and trim and pad
    analysis_arg_parsing.analysis_group_parsing(args.reference_fasta,args.trim_start,args.trim_end,args.catchment_size,args.downsample,args.query_limit,config)

    # Checks background data exists and is the right format.
    # Checks same number of records supplied for csv, fasta and (optional) SNP file. 
    data_arg_parsing.data_group_parsing(args.debug,args.datadir,args.background_csv,args.background_SNPs,args.background_fasta,args.background_tree,args.background_column,args.fasta_column,config)

    # Sort out where the query info is coming from, csv or id string, optional fasta seqs.
    # Checks if they're real files, of the right format and that QC args sensible values.
    input_arg_parsing.input_query_parsing(args.input_csv,args.input_column,args.ids,args.from_metadata,config)
    input_arg_parsing.input_fasta_parsing(args.fasta,args.max_ambiguity,args.min_length,config)

    input_arg_parsing.from_metadata_parsing(config)
    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    snakefile = data_install_checks.get_snakefile(thisdir)


    # Checks there are records to run and matches them up from background or supplied fasta
    # merges the metadata to a master metadata
    # runs supplied fasta qc
    query_metadata, passed_qc_fasta, found_in_background_data = input_data_parsing.query_check_against_background_merge_input(config)

    # Define what's going to go in the report and sort global report options 
    # stored under config = { "report_content": [1, 2, 3, 4], "reports": [1,2,3,4],[1,2]}
    name_dict = report_arg_parsing.parse_global_report_options(args.report_content, args.report_column, args.anonymise, args.date_column, args.background_date_column, args.location_column, config)
    report_arg_parsing.parse_optional_report_content(args.table_content, args.timeline_dates, args.timeline_colours, config)
    report_arg_parsing.parse_map_options(args.background_map_date_restriction, args.background_map_location, args.map_file, args.longitude_column, args.latitude_column, found_in_background_data, config)

    # sets up the output dir, temp dir, and data output desination
    directory_setup.output_group_parsing(args.outdir, args.output_prefix, args.overwrite,args.datestamp, args.output_data, args.tempdir, args.no_temp, config)

    # write the merged metadata, the extracted passed qc supplied fasta and the extracted matched fasta from the background data
    input_data_parsing.write_parsed_query_files(query_metadata,passed_qc_fasta,found_in_background_data, config)

    if config["anonymise"]:
        global_report_functions.write_anon_names_to_file(config, name_dict)

    # ready to run? either verbose snakemake or quiet mode

    if config["verbose"]:
        print(red("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config["tempdir"],config=config, cores=config["threads"],lock=False
                                    )
    else:
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=config["tempdir"],
                                    config=config, cores=config["threads"],lock=False,quiet=True,log_handler=config["log_api"]
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1
    

if __name__ == '__main__':
    main()