#!/usr/bin/env python3

import argparse
from civet.utils.log_colours import green,cyan,red
from civet.report_functions import report
import os
import yaml
from Bio import SeqIO

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Render report from yaml file.')

    parser.add_argument("--csv", action="store", type=str, dest="csv")
    parser.add_argument("--report", action="store",dest="report")
    parser.add_argument("--yaml", action="store", type=str, dest="yaml")
    return parser.parse_args()

def render_report():
    args = parse_args()

    with open(args.yaml, 'r') as f:
        config_loaded = yaml.safe_load(f)
    
    
    for report_to_generate in args.report.split(","):
            report.make_report(args.csv,report_to_generate,config_loaded)

if __name__ == '__main__':

    render_report()