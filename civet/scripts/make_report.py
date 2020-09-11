#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys
from input_qc_functions import parse_yaml_file
import yaml
thisdir = os.path.abspath(os.path.dirname(__file__))

def parse_args():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("--config", required=True, help="config yaml file", dest="config") 
    return parser.parse_args()

def make_report():
    args = parse_args()

    config = {}

    parse_yaml_file(args.config,config)

    with open(config["outfile"], 'w') as pmd_file:
        change_line_dict = {}

        for key in config:
            change_line_dict[key] = f'{key} = "{config[key]}"\n'
        
        with open(config["report_template"]) as f:
            for l in f:
                if "##CHANGE" in l:
                    for key in change_line_dict:
                        if key in l:
                            new_l = change_line_dict[key]
                else:
                    new_l = l

                pmd_file.write(new_l)
    
    weave(config["outfile"], doctype = "pandoc", figdir=config["rel_figdir"])
    
if __name__ == "__main__":
    make_report()
