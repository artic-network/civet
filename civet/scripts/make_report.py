#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys
from reportfunk.funks import io_functions as qcfunk

import yaml
thisdir = os.path.abspath(os.path.dirname(__file__))

def parse_args():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("--config", required=True, help="config yaml file", dest="config") 
    return parser.parse_args()

def make_free_text_dict(config):

    free_text_dict = {}

    free_text_dict["##INSERT_TITLE"] = config["title"] + "\n"
    free_text_dict["##OUTBREAKID"] = config["outbreak_id"] + "\n"
    free_text_dict["##DATE"] = config["report_date"] + "\n"
    free_text_dict["##AUTHORS"] = config["authors"] + "\n"
    free_text_dict["##DESCRIPTION"] = config["description"]
    free_text_dict["##CONCLUSIONS"] = config["conclusions"]


    config.pop('title', None)
    config.pop('report_date', None)
    config.pop('authors', None)
    config.pop('description', None)
    config.pop('conclusions', None)

    return free_text_dict

def make_report():
    
    args = parse_args()
    config = {}
    qcfunk.parse_yaml_file(args.config,config)
    free_text_dict = make_free_text_dict(config)

    path = os.path.join(config["outdir"], "report", "figures","svg_figures")
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    
    config["svg_figdir"] = path

    with open(config["outfile"], 'w') as pmd_file:
        change_line_dict = {}

        for key, value in config.items():
            if type(value) == bool:
                new_value = f'{key} = {value}\n'
            else:
                new_value = f'{key} = "{value}"\n'
            
            change_line_dict[key] = new_value

        with open(config["report_template"]) as f:
            for l in f:
                line_written = False
                if "##INSERT_ARGUMENTS" in l:
                        pmd_file.write("".join(list(change_line_dict.values())))
                        line_written = True

                elif "##APPENDIX" in l:
                    if config["omit_appendix"]:
                        new_l = ""
                        line_written = True
                    else:
                        line_written = True
                        with open(config['appendix']) as f:
                            for l in f:
                                pmd_file.write(l)
                else:
                    for k,v in free_text_dict.items():
                        if k in l:
                            if "'''" not in v:
                                new_l = str(v)
                            else:
                                new_l = v.lstrip("'").rstrip("'")

                            pmd_file.write(new_l)
                            line_written = True

                if not line_written:
                    pmd_file.write(l)
        
    
    weave(config["outfile"], doctype = "pandoc", figdir=config["figdir"])
    
if __name__ == "__main__":
    make_report()
