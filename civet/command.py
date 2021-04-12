#!/usr/bin/env python3
from civet import __version__

from civet.input_parsing import initialising as init
from civet.input_parsing import input_arg_parsing
from civet.input_parsing import data_arg_parsing

from civet.utils import misc

import os
import sys
import argparse

cwd = os.getcwd()

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
    i_group.add_argument('-l','--min-length', action="store", type=int,help="Minimum query length allowed to attempt analysis. Default: 10000",dest="min_length")
    
    i_group.add_argument('-fm','--from-metadata',nargs='*', dest="from_metadata",help="Generate a query from the metadata file supplied. Define a search that will be used to pull out sequences of interest from the large phylogeny. E.g. -fm adm2=Edinburgh sample_date=2020-03-01:2020-04-01")

    d_group = parser.add_argument_group('Background data options')
    d_group.add_argument('-d','--datadir', action="store",help="Directory containing the background data files.")
    d_group.add_argument("-bc","--background-csv",action="store",dest="background_csv",help="Custom metadata file that corresponds to the large global tree/ alignment. Should have a column `sequence_name`.")
    d_group.add_argument("-bf","--background-fasta", action="store", dest="background_fasta", help="Custom background fasta file.")
    d_group.add_argument("-dcol",'--data-column', action="store",help="Column in background csv to match with input IDs. Default: sequence_name", dest="data_column")

    o_group = parser.add_argument_group('Output options')
    o_group.add_argument('-o','--output-prefix',action="store",help="Prefix of output directory & report name: Default: civet",dest="output_prefix")
    o_group.add_argument('--outdir', action="store",help="Output directory. Default: current working directory")
    
    o_group.add_argument('--output-data',action="store_true",help="Output intermediate data files to the output directory",dest="output_data")
    
    o_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files. For development/ debugging purposes.",dest="no_temp")
    o_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")

    
    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument("--safety-level", action="store", type=int, dest="safety_level",help="Level of anonymisation for users. Options: 0 (no anonymity), 1 (no COGIDs on background data), 2 (no adm2 on data). Default: 1")
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


    config = init.setup_config_dict(cwd, args.config)


    input_arg_parsing.input_query_parsing(args.input_csv,args.input_column,args.ids,config)


    input_arg_parsing.input_fasta_parsing(args.fasta,args.max_ambiguity,args.min_length,config)

    data_arg_parsing.data_group_parsing(args.datadir,args.background_csv,args.background_fasta,args.data_column,config)

    for i in sorted(config):
        print(i, config[i])

if __name__ == '__main__':
    main()