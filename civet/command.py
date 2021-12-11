#!/usr/bin/env python3
from civet import __version__

from civet.input_parsing import initialising as init
from civet.input_parsing import input_arg_parsing
from civet.input_parsing import data_arg_parsing
from civet.input_parsing import analysis_arg_parsing
from civet.input_parsing import input_data_parsing
from civet.input_parsing import report_arg_parsing
from civet.input_parsing import generate_background_parsing

from civet.output_options import directory_setup
from civet.report_functions import global_report_functions

from civet.utils import misc
from civet.utils import dependency_checks
from civet.utils import data_install_checks

import civet.utils.custom_logger as custom_logger
from civet.utils.log_colours import green,cyan,red
from civet.utils.config import *

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
    i_group.add_argument('-ids',"--id-string", action="store",help="Comma-separated id string with one or more query ids. Example: `EDB3588,EDB3589`", dest="id_string")
    i_group.add_argument('-i',"--input-metadata", action="store",help="Input csv file with minimally an <input_id_column> header. Default=`name`", dest="input_metadata")
    i_group.add_argument('-f','--input-sequences', action="store",help="Optional fasta file. Sequence IDs must match to a query ID specified either in the input csv or ID string", dest="input_sequences")
    i_group.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the background data. Example: -fm country=Ireland sample_date=2020-03-01:2020-04-01")
    i_group.add_argument('-mq','--max-queries', type=int, action="store",dest="max_queries",help="Max number of queries. Default: `5000`")
    i_group.add_argument('--focal_alignment', action="store",
                            help="Optional alignment of focal sequences for global snipit", dest="focal_alignment")
    i_group.add_argument('--reference_name', action="store",
                         help="Optional input for the reference name in the focal alignment. Default: Reference",
                         dest="reference_name")

    ic_group = parser.add_argument_group('Input column configuration')
    ic_group.add_argument('-icol',"--input-id-column", action="store", dest="input_id_column",help="Column in input csv file to match with database. Default: `name`")
    ic_group.add_argument("-idisp", "--input-display-column", action="store", dest="input_display_column", help="Column containing alternative names to display in the report. Example: `patient_IDs`")
    ic_group.add_argument("-idate","--input-date-column", action="store", dest="input_date_column", help="Column in input query with date data in. Default: `sample_date`")

    is_group = parser.add_argument_group("Input sequences options")
    is_group.add_argument('-n','--max-ambiguity', action="store", type=float,help="Maximum proportion of Ns allowed to attempt analysis. Default: `0.5`",dest="max_ambiguity")
    is_group.add_argument('-l','--min-length', action="store", type=int,help="Minimum query length allowed to attempt analysis. Default: `20000`",dest="min_length")
    is_group.add_argument('-ts','--trim-start',action="store", type=int, dest="trim_start",help="Genome position to trim and pad to when aligning input sequences. Default: `265`")
    is_group.add_argument('-te','--trim-end', action="store",type=int, dest="trim_end",help="Genome position to trim and pad from when aligning input sequences. Default: `29674`")

    d_group = parser.add_argument_group('Background data options')
    d_group.add_argument('-d','--datadir', action="store",help="Directory containing the background data files.")
    d_group.add_argument("-bm","--background-metadata",action="store",dest="background_metadata",help="Custom metadata file for all background data. Should have a column matching <-bicol/--background-id-column>")
    d_group.add_argument("-bsnp","--background-snps",action="store",dest="background_snps",help="Optional SNP file for all background data. Civet will calculate this file if not supplied, which may take some time")
    d_group.add_argument("-bseq","--background-sequences", action="store", dest="background_sequences", help="Custom background sequence file for all background data. Sequence IDs should match the background metadata id column, or specify another column using <-biseq/--background-sequence-id>")
    d_group.add_argument("-bt","--background-tree", action="store", dest="background_tree", help="Custom background tree file for all background data. Tip names should match the background metadata background_column. *Coming soon*")
    
    bc_group = parser.add_argument_group('Background column configuration')
    bc_group.add_argument("-bicol",'--background-id-column', action="store", dest="background_id_column",help="Column in background metadata to match with input IDs. Default: `sequence_name`")
    bc_group.add_argument("-sicol",'--sequence-id-column', action="store", dest="sequence_id_column",help="Column in background data to match with sequence IDs. Default: <-bicol/--background-id-column>")
    bc_group.add_argument("-bdate", "--background-date-column", action="store", dest="background_date_column", help="Column with date information in background metadata. Default: `sample_date`")
    bc_group.add_argument("-bloc","--background-location-column", action="store",dest="background_location_column", help="Column containing geographical data in the background metadata. Default: `country`")

    o_group = parser.add_argument_group('Output options')
    o_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: `civet-2021-XX-YY`")
    o_group.add_argument('-p','--output-prefix',action="store",help="Prefix of output directory & report name: Default: `civet`",dest="output_prefix")
    o_group.add_argument('--datestamp', action="store",help="Append datestamp to directory name when using <-o/--outdir>. Default: <-o/--outdir> without a datestamp")
    o_group.add_argument('--overwrite', action="store_true",help="Overwrite output directory. Default: append an incrementing number if <-o/--outdir> already exists")
    o_group.add_argument('--output-data',action="store_true",help="Output intermediate data files to the output directory",dest="output_data")
    o_group.add_argument('-temp','--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: `$TMPDIR`")
    o_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files. For development/ debugging purposes",dest="no_temp")

    dc_group = parser.add_argument_group('Background data curation')
    dc_group.add_argument("-bd","--generate-civet-background-data",dest="generate_civet_background_data",action="store",help="A sequence file to create background metadata, alignment and SNP file from.")
    dc_group.add_argument("--background-data-checks",dest="debug",action="store_true",help="Run checks on custom background data files, not run by default")
    dc_group.add_argument("--background-data-outdir",dest="background_data_outdir",action="store",help="Directory to output the civet background data. Default: `civet_data`")
    dc_group.add_argument("--primary-field-delimiter",dest="primary_field_delimiter",action="store",help="Primary sequence header field delimiter to create metadata file from. Default: `|`")
    dc_group.add_argument("--primary-metadata-fields",dest="primary_metadata_fields",action="store",help="Primary sequence header fields to create metadata file from. Default: `sequence_name,gisaid_id,sample_date`")
    dc_group.add_argument("--secondary-fields",dest="secondary_fields",action="store_true",help="Parse header for secondary sequence fields to add to metadata file. Default: `False`")
    dc_group.add_argument("--secondary-field-delimiter",dest="secondary_field_delimiter",action="store",help="Secondary sequence header field delimiter to create metadata file from. Default: `/`")
    dc_group.add_argument("--secondary-field-location",dest="secondary_field_location",action="store",help="Secondary sequence header location within primary field list. Default: `0` (i.e. the first field)")
    dc_group.add_argument("--secondary-metadata-fields",dest="secondary_metadata_fields",action="store",help="Secondary sequence header fields to create metadata file from. Default: `virus,country,sequence_id,year`")

    c_group = parser.add_argument_group("Catchment options")
    c_group.add_argument('-snpd','--snp-distance', type=int, action="store",dest="snp_distance",help="Define radius of catchment by number of SNPs from query. Default: `2`")
    c_group.add_argument('-pushd','--push-distance', type=int, action="store",dest="push_distance",help="Define node radius to push from query. Default: `2`")
    c_group.add_argument('--snp-distance-up', type=int, action="store",dest="snp_distance_up",help="Define radius of catchment by number of SNPs from query. Default: <-snpd/--snp-distance>")
    c_group.add_argument('--snp-distance-down', type=int, action="store",dest="snp_distance_down",help="Define radius of catchment by number of SNPs from query. Default: <-snpd/--snp-distance>")
    c_group.add_argument('--snp-distance-side', type=int, action="store",dest="snp_distance_side",help="Define radius of catchment by number of SNPs from query. Default: <-snpd/--snp-distance>")
    c_group.add_argument('-cs','--catchment-background-size', type=int, action="store",dest="catchment_background_size",help="Max number of background sequences in a catchment. Catchments with more background sequences than this will be downsampled prior to tree building. Default: `100`")
    c_group.add_argument('-ds','--downsample', nargs='*', action="store",dest="downsample",help="""Configuration of catchment downsampling. Options: `random`, `enrich` or `normalise`. Default: `random`).
If using enrich mode, indicate the factor (Default: `10`), the column name and the field to enrich for.   
Example: --downsample mode=enrich factor=10 sample_date=2021-02-04:2021-03-04   
If using normalise mode, indicate the column to normalise across.   
Example: --downsample mode=normalise country""")

    t_group = parser.add_argument_group("Tree options")
    t_group.add_argument("-ta","--tree-annotations", action="store", dest="tree_annotations", help="Comma separated string of metadata columns to annotate catchment trees with, can then be displayed in the report.")
    t_group.add_argument("-mt","--max-tree-size", action="store", dest="max_tree_size", help="Maximum number of sequences allowed in a given catchment tree. Catchments with more than this will be summarised in a table rather than a tree. Default: 500")

    r_group = parser.add_argument_group("Report configuration options")
    r_group.add_argument("-rp", "--report-preset", nargs='*', action="store", dest="report_preset", help="""Specify one or more preset options to configure report content, full report configuration options available using <-rc/--report-content>   
Preset options: `the_usual` (1,2,3,4,5), `the_works` (1,2,3,4,5,7), `the_whole_shebang` (1,2,3,4,5,6,7), `hold_the_sauce` (1,2)   
Default: `the_usual`""")
    r_group.add_argument("-rc", "--report-content", nargs='*', action="store", dest="report_content", help="""One or more comma separated numeric strings to define the report content. Default: `1,2,3,4,5`    
1: Query summary tables, 
2: Catchment summary tables, 
3: Catchment trees, 
4: snipit plots, 
5: Timeline, 
6: Background map, 
7: Query map,
8: Time series""")
    r_group.add_argument('--global_snipit', action="store_true",
                              help="Create a global snipit figure for all focal sequences", dest="global_snipit")
    r_group.add_argument("-rt", "--report-title", action="store", dest="report_title", help="""Title to display in report. Default: `civet report`""")
    r_group.add_argument("--anonymise", action="store_true", dest="anonymise",help="Generates arbitrary labels for sequences for dissemination")
    r_group.add_argument("-ct", "--colour-theme", action="store", dest="colour_theme", help="""Report theme colour. Default: `#7178bc`""")
    r_group.add_argument("-cmap", "--colour-map", action="store", dest="colour_map", help="""Comma separated string of hex codes or names of HTML compatible colours to colour factors in report by. Default: `"#86b0a6,#565b7b,#e9c46a,#e1998a,#d19c2c,#264653,#4888a3,#c77939,#b56576,#eb7e83,#F46D43,#9e2a2b,#84a98c"`""")   

    tb_group = parser.add_argument_group("Table options (report options 1 and 2)")
    tb_group.add_argument("--query-table-content", action='store', dest="query_table_content", help="Columns to include in the table for queries. Default: `<-bicol/--background-id_column>,<-bdate/--background_date_column>,source,lineage,country,catchment`")
    tb_group.add_argument("--catchment-table-content", action='store', dest="catchment_table_content", help="Columns to include in the summary table for catchments. Default: `count,country,lineage`")
    tb_group.add_argument("--mutations",action="store",dest="mutations",help="Comma separated string of mutations to type query sequences and display in query summary table.")

    tl_group = parser.add_argument_group("Timeline options (report option 5)")
    tl_group.add_argument("-tdate", "--timeline-dates", action='store', dest="timeline_dates", help="Data to generate a timeline as a comma separated string. Default: <-idate/--input-date-column>")
    tl_group.add_argument("-tgc", "--timeline-group-column", action='store', dest="timeline_group_column", help="Column to group sequences by to show as single lines on the timeline figures. Default: input_display_column")

    bm_group = parser.add_argument_group("Background map options (report option 6)")
    bm_group.add_argument("-bmfile","--background-map-file", action="store", dest="background_map_file", help="Topojson containing polygons to plot background diversity on. Must be an online resource eg on a Github pages website.")
    bm_group.add_argument("--centroid-file", action="store", dest="centroid_file", help="csv containing centroids matching locations in background_map_file. Must be provided if custom geojson/json is provided for background mapping. Headers must be location, latitude and longitude.")
    bm_group.add_argument("-bmloc,", "--background-map-location", action="store", dest="background_map_location", help="Comma separated list containing locations to show background lineage diversity for. Default is all locations at the appropriate administrative level.")
    bm_group.add_argument("-bmdr","--background-map-date-range", dest="background_map_date_range", action="store", help="Date range for mapping background lineage diversity. Can be an integer (number of days either side of queries to restrict to) or a date range, format='YYYY-MM-DD:YYYY-MM-DD'")
    bm_group.add_argument("-bmcol","--background-map-column", dest="background_map_column", action="store", help="Column in background metadata containing location to map background lineage diversity by.")
    bm_group.add_argument("--background-map-colours", dest="background_map_colours", action="store", help="list of 20 colours to colour most common lineages by in background map.")
    bm_group.add_argument("--background-map-other-colours", dest="background_map_other_colours", action="store", help="Comma separated string of two colours to colour lineages that make up less than 5%% of a country's sequences, and those that no tin the 20 most common lineages in the whole dataset.")

    qm_group = parser.add_argument_group("Query map options (report option 7)")
    qm_group.add_argument("-qmfile","--query-map-file", action="store", dest="query_map_file", help="Topojson containing polygons to plot queries on. Must be an online resource eg on a Github pages website.")
    qm_group.add_argument("-lat","--latitude-column", dest="latitude_column", action="store", help="Column containing latitude coordinate information to plot queries on a map")
    qm_group.add_argument("-long","--longitude-column", dest="longitude_column", action="store", help="Column containing longitude coordinate information to plot queries on a map")

    bs_group = parser.add_argument_group("Query time series options (report option 8)")
    bs_group.add_argument("--series-colour-factor", action='store', dest="series_colour_factor", help="Comma separated string of columns to colour time series by. Default: lineage")

    misc_group = parser.add_argument_group('Misc options')
    misc_group.add_argument("--civet-mode", action="store", dest='civet_mode', help="If CLIMB then import UK specific modules. Default=`GLOBAL`")
    misc_group.add_argument("-r","--reference-sequence",action="store",dest="reference_sequence",help="Custom reference genome to map and pad against. Must match the reference the background sequence alignment was generated from.")
    misc_group.add_argument("--art",action="store_true",help="Print art")
    misc_group.add_argument("--acknowledgements",action="store_true",help="Print thanks")
    misc_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
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
        misc.be_arty()
        sys.exit(0)

    if args.acknowledgements:
        misc.full_acknowledgements()
        sys.exit(0)

    dependency_checks.check_dependencies(dependency_list, module_list)
    if args.mutations:
        dependency_checks.check_scorpio_mutations(args.mutations)
    
    # Initialise config dict
    config = init.setup_config_dict(cwd,args.config)

    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    snakefile = data_install_checks.get_snakefile(thisdir)
    
    # Threads and verbosity to config
    init.misc_args_to_config(args.verbose,args.threads,args.civet_mode,config)
    init.set_up_verbosity(config)

    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    
    if args.generate_civet_background_data:
        generate_background_parsing.parse_generate_background_args(args.generate_civet_background_data,args.background_data_outdir,args.primary_field_delimiter,args.primary_metadata_fields,args.secondary_fields,args.secondary_field_delimiter,args.secondary_field_location,args.secondary_metadata_fields,config)
        snakefile = data_install_checks.get_generator_snakefile(thisdir)
        if config[KEY_VERBOSE]:
            print(red("\n**** CONFIG ****"))
            for k in sorted(config):
                print(green(f" - {k}: ") + f"{config[k]}")
            status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        config=config, cores=config[KEY_THREADS],lock=False
                                        )
        else:
            status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,
                                        config=config, cores=config[KEY_THREADS],lock=False,quiet=True,log_handler=config[KEY_LOG_API]
                                        )
        if status: # translate "success" into shell exit code of 0
            return 0   

        return 1


    # Checks background data exists and is the right format.
    # Checks same number of records supplied for csv, fasta and (optional) SNP file. 
    data_arg_parsing.data_group_parsing(args.debug,args.datadir,args.background_metadata,args.background_snps,args.background_sequences,args.background_tree,args.background_id_column,args.sequence_id_column,config)

    # Analysis options, including ref and trim and pad
    analysis_arg_parsing.analysis_group_parsing(args.reference_sequence,args.trim_start,args.trim_end,args.max_queries,config)

    analysis_arg_parsing.catchment_group_parsing(args.catchment_background_size,args.downsample,args.push_distance,args.snp_distance,args.snp_distance_up,args.snp_distance_down,args.snp_distance_side,config)

    # Sort out where the query info is coming from, csv or id string, optional fasta seqs.
    # Checks if they're real files, of the right format and that QC args sensible values.
    input_arg_parsing.input_query_parsing(args.input_metadata,args.input_id_column,args.id_string,args.from_metadata,config)
    input_arg_parsing.input_fasta_parsing(args.input_sequences,args.max_ambiguity,args.min_length,config)

    input_arg_parsing.from_metadata_parsing(config)

    snakefile = data_install_checks.get_snakefile(thisdir)

    # Checks there are records to run and matches them up from background or supplied fasta
    # merges the metadata to a master metadata
    # runs supplied fasta qc
    query_metadata, passed_qc_fasta, found_in_background_data = input_data_parsing.query_check_against_background_merge_input(config)

    # Define what's going to go in the report and sort global report options 
    # stored under config = { "report_content": [1, 2, 3, 4], "reports": [1,2,3,4],[1,2]}
    name_dict = report_arg_parsing.parse_global_report_options(args.report_content,args.report_preset, args.input_display_column, args.anonymise, args.input_date_column, args.background_date_column, args.background_location_column, config)
    report_arg_parsing.parse_optional_report_content(args.query_table_content,args.mutations, args.timeline_dates, args.timeline_group_column, args.colour_theme, args.colour_map, config)
    report_arg_parsing.parse_map_options(args.background_map_date_range, args.background_map_column, args.background_map_file, args.centroid_file, args.background_map_location, args.query_map_file, args.longitude_column, args.latitude_column, found_in_background_data, args.background_map_colours, args.background_map_other_colours,config)
    report_arg_parsing.parse_tree_options(args.tree_annotations,args.max_tree_size, config)
    report_arg_parsing.parse_series_options(args.series_colour_factor,query_metadata, config[KEY_INPUT_DATE_COLUMN],config)

    # sets up the output dir, temp dir, and data output desination
    directory_setup.output_group_parsing(args.outdir, args.output_prefix, args.overwrite, args.datestamp, args.output_data, args.tempdir, args.no_temp, config)

    # write the merged metadata, the extracted passed qc supplied fasta and the extracted matched fasta from the background data
    input_data_parsing.write_parsed_query_files(query_metadata,passed_qc_fasta,found_in_background_data, config)

    if config[KEY_ANONYMISE]:
        global_report_functions.write_anon_names_to_file(config, name_dict)

    # ready to run? either verbose snakemake or quiet mode
    for k in sorted(config):
        print(green(f" - {k}: ") + f"{config[k]}")

    if config[KEY_VERBOSE]:
        print(red("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR],config=config, cores=config[KEY_THREADS],lock=False
                                    )
    else:
        
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=config[KEY_TEMPDIR],
                                    config=config, cores=config[KEY_THREADS],lock=False,quiet=True,log_handler=config[KEY_LOG_API]
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1



if __name__ == '__main__':
    main()
