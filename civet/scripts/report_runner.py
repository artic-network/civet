#!/usr/bin/env python3

import argparse
from civet.utils.log_colours import green,cyan,red
from civet.report_functions import report
import os
import yaml
from Bio import SeqIO
import sys

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Render report from yaml file.')

    parser.add_argument("--csv", action="store",help="Metadata file to generate report with.", type=str, dest="csv")
    parser.add_argument("--report", action="store", help="Name of report to generate.",dest="report")
    parser.add_argument("--configfile", action="store", type=str,help="Yaml file to generate report with.", dest="configfile")
    parser.add_argument("--config", action="store", nargs='*',help="Overwrite any options in the config yaml file.",,dest="config")
    return parser.parse_args()

def render_report():
    args = parse_args()
    
    with open(args.configfile, 'r') as f:
        config_loaded = yaml.safe_load(f)
    
    if args.config:
        for arg in args.config:
            if '=' not in arg:
                sys.stderr.write(cyan(f'Error: please enter `--config` items as key=value pairs.\n') + f'Invalid: {arg}\n')
                sys.exit(-1)
            else:
                key,value = arg.split("=")
                config_loaded[key] = value

    for report_to_generate in args.report.split(","):
            report.make_report(args.csv,report_to_generate,config_loaded)

if __name__ == '__main__':

    render_report()