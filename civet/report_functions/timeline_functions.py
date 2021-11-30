import csv
from collections import defaultdict
import json
import os
import sys
from civet.utils.log_colours import cyan 
from civet.utils import misc

from civet.report_functions import global_report_functions

def timeline_checking(timeline_dates, timeline_group_column, config):  

    """
    parses the report group argument
    --timeline-dates
    --timeline_group_column
    """

    misc.add_arg_to_config("timeline_dates",timeline_dates,config) #default is None
    misc.add_arg_to_config("timeline_group_column", timeline_group_column, config)

    if config["timeline_dates"] or config["timeline_group_column"]:

        if "input_metadata" in config:
            with open(config["input_metadata"]) as f:
                data = csv.DictReader(f)
                input_headers = data.fieldnames
        else:    
            input_headers = []

        with open(config["background_metadata"]) as f:
            data = csv.DictReader(f)
            background_headers = data.fieldnames

    if config["timeline_group_column"]:
        if config["timeline_group_column"] not in input_headers and config["timeline_group_column"] not in background_headers:
            sys.stderr.write(cyan(f"Error: {config['timeline_group_column']} (specified for use for grouping the timeline plot) not found in either metadata file.\n"))
            sys.exit(-1)
    else:
        config["timeline_group_column"] = config["input_display_column"]

    if config["timeline_dates"]:

        date_header_list = config["timeline_dates"].split(",")

        input_cols_to_check = []
        background_cols_to_check = []

        for header in date_header_list:
            if header not in input_headers and header not in background_headers:
                sys.stderr.write(cyan(f"Error: {header} (specified for use in timeline plot) not found in either metadata file.\n"))
                sys.exit(-1)

            elif header in input_headers:
                input_cols_to_check.append(header)
            elif header in background_headers:
                background_cols_to_check.append(header)

        count = 0
        if "input_metadata" in config:
            with open(config["input_metadata"]) as f:
                data = csv.DictReader(f)
                for l in data:
                    for header in input_cols_to_check:
                        count += 1
                        if l[header] != "":
                            misc.check_date_format(l[header], count, header)
        
        count = 0
        with open(config["background_metadata"]) as f:
            data = csv.DictReader(f)
            for l in data:
                for header in background_cols_to_check:
                    count += 1
                    if l[header] != "":
                        misc.check_date_format(l[header], count, header)


    elif config["input_date_column"] and config["background_date_column"]:
        if config["input_date_column"] != config["background_date_column"]:
            config["timeline_dates"] = f'{config["input_date_column"]},{config["background_date_column"]}'
        else:
            config["timeline_dates"] = config["input_date_column"]
    elif not config["background_date_column"] and config["input_date_column"]:
        config["timeline_dates"] = config["input_date_column"]
    elif config["background_date_column"] and not config["input_date_column"]:
        config["timeline_dates"] = config["background_date_column"]
    
    else:
        sys.stderr.write(cyan(f"Error: Timeline option given in report content argument, but no data provided. Please use -td/--timeline-dates to specify column headers containing data.\n"))
        sys.exit(-1)


    return config


def make_timeline_json(catchment,config):

    if type(config["timeline_dates"]) == str:
        date_cols = config["timeline_dates"].split(",")
        config["timeline_dates"] = date_cols

    date_cols = config["timeline_dates"]

    overall = defaultdict(dict)
    catchment_string = f'{catchment}_timeline'
    dict_list = []
    with open(config["query_metadata"]) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['catchment'] == catchment:
                for col in date_cols:
                    new_dict = {}
                    new_dict["sequence_name"] = l[config["timeline_group_column"]]
                    new_dict["date"] = l[col]
                    new_dict["date_type"] = col
                    dict_list.append(new_dict)
                
    overall[catchment_string] = dict_list  

    return json.dumps(overall)
