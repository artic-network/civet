#!/usr/bin/env python3
from civet.utils import log_colours as colour
import sys
import os

def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

## i_group
def i_group_parsing(input_csv,ids,config):

    add_arg_to_config("ids",ids,config)

    add_arg_to_config("input_csv",input_csv,config)

    if "ids" in config and "input_csv" in config:

        sys.stderr.write(cyan(f"Error: it looks like you've provide a csv file and an ID string, please provide one or the other.\n"))
        sys.exit(-1)
    elif "ids" in config:
        if not type(config["ids"]) is list:
            config["ids"] = config["ids"].split(",")
        
        config["ids"] = list(set(config["ids"]))

    print("Unique IDs input:", len(config["ids"]))
