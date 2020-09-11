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

# def generate_title - so if they provide sc but no title, then make it here
# def make_free_text_dict

def make_report(inputs, report_args_file):
    
    arguments = get_report_arguments(report_args_file)

    config = {}

    parse_yaml_file(args.config,config)

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
        
        #make dict with these, add in \ns on the end
        title = "# testy testy title\n"
        date = "2020-08-01\n"
        description = "We're looking at xyz investigation\n"
        authors = "- Martin \n - Matt \n"


        #make_pmd_file(pmd_file, md_template, change_line_dict, free_text)
        free_text_dict = {}

        #will have to work out how to generate this from the incoming config file
        free_text_dict["##INSERT_TITLE"] = title 
        free_text_dict["##DATE"] = f'This investigation was started on {date}'
        free_text_dict["##DESCRIPTION"] = description
        free_text_dict["##AUTHORS"] = authors


        
        with open(config["report_template"]) as f:
            for l in f:
                line_written = False
                if "##INSERT_ARGUMENTS" in l:
                        pmd_file.write("".join(list(change_line_dict.values())))
                        line_written = True
                else:
                    for k,v in free_text_dict.items():
                        if k in l:
                            new_l = str(v)
                            pmd_file.write(new_l)
                            line_written = True

                if not line_written:
                    pmd_file.write(l)

       
        
    
    weave(config["outfile"], doctype = "pandoc", figdir=config["figdir"])
    
if __name__ == "__main__":
    make_report()
