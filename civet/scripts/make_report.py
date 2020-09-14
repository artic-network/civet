#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys
from input_qc_functions import parse_yaml_file
import yaml
thisdir = os.path.abspath(os.path.dirname(__file__))

def get_report_arguments(arg_file):

    args = []

    with open(arg_file) as f:
        for l in f:
            arg = l.strip("\n")
            arg = arg.replace("-", "_")
            args.append(arg)

    return args

def parse_args():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("--config", required=True, help="config yaml file", dest="config") 
    return parser.parse_args()

def make_free_text_dict(config):

    free_text_dict = {}

    free_text_dict["##INSERT_TITLE"] = config["title"] + "\n"
    free_text_dict["##DATE"] = config["report_date"] + "\n"
    free_text_dict["##AUTHORS"] = config["authors"] + "\n"
    free_text_dict["##DESCRIPTION"] = config["description"]
    
    config.pop('title', None)
    config.pop('report_date', None)
    config.pop('authors', None)
    config.pop('description', None)

    return free_text_dict

def make_report():
    
    args = parse_args()
    
    config = {}
    config = parse_yaml_file(args.config,config)
    free_text_dict = make_free_text_dict(config)

    with open(config["outfile"], 'w') as pmd_file:
        change_line_dict = {}

        for key, value in config.items():
            if key == "cog_metadata":
                new_key = "full_metadata_file"
            else:
                new_key = key

            new_value = f'{new_key} = "{value}"\n'
            change_line_dict[new_key] = new_value
        
        if config["add_bars"]:
            change_line_dict["add_bars"] = f'add_bars = "{config["add_bars"]}"\n'
        else:
            change_line_dict["add_bars"] = 'add_bars = ""\n'

        with open(config["report_template"]) as f:
            for l in f:
                line_written = False
                if "##INSERT_ARGUMENTS" in l:
                        pmd_file.write("".join(list(change_line_dict.values())))
                        line_written = True
                else:
                    for k,v in free_text_dict.items():
                        if k in l:
                            if k != "##DESCRIPTION":
                                new_l = str(v)
                            else:
                                new_l = v

                            pmd_file.write(new_l)
                            line_written = True

                if not line_written:
                    pmd_file.write(l)
        
    
    weave(config["outfile"], doctype = "pandoc", figdir=config["figdir"])
    
if __name__ == "__main__":
    make_report()
