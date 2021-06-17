import csv
from collections import defaultdict
import json
import os
import sys
from civet.utils.log_colours import cyan 
from civet.utils import misc


def timeline_checking(timeline_dates, config):  

    """
    parses the report group argument
    --timeline-dates
    """

    misc.add_arg_to_config("timeline_dates",timeline_dates,config) #default is None

    if config["timeline_dates"]:

        date_header_list = config["timeline_dates"].split(",")

        if "input_csv" in config:
            with open(config["input_csv"]) as f:
                data = csv.DictReader(f)
                input_headers = data.fieldnames
        else:    
            input_headers = []

        with open(config["background_csv"]) as f:
            data = csv.DictReader(f)
            background_headers = data.fieldnames

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
        if "input_csv" in config:
            with open(config["input_csv"]) as f:
                data = csv.DictReader(f)
                for l in data:
                    for header in input_cols_to_check:
                        count += 1
                        if l[header] != "":
                            misc.check_date_format(l[header], count, header)
        
        count = 0
        with open(config["background_csv"]) as f:
            data = csv.DictReader(f)
            for l in data:
                for header in background_cols_to_check:
                    count += 1
                    if l[header] != "":
                        misc.check_date_format(l[header], count, header)


    else:
        sys.stderr.write(cyan(f"Error: Timeline option given in report content argument, but no data provided. Please use -td/--timeline-dates to specify column headers containing data.\n"))
        sys.exit(-1)


    return config


def make_timeline_json(config):

    date_cols = config["timeline_dates"].split(",")
    config["timeline_dates"] = date_cols

    overall = defaultdict(dict)

    with open(config["query_metadata"]) as f:
        data = csv.DictReader(f)
        for l in data:

            catchment_string = f"{l['catchment']}_timeline"
           
            if catchment_string in overall:
                dict_list = overall[catchment_string]
            else:
                dict_list = []
            
            for col in date_cols:
                new_dict = {}
                new_dict["sequence_name"] = l[config["report_column"]]
                new_dict["date"] = l[col]
                new_dict["date_type"] = col
                dict_list.append(new_dict)
            
            overall[catchment_string] = dict_list  

    timeline_json = os.path.join(config["tempdir"],'timeline_data.json')

    with open(timeline_json, 'w') as outfile:
        json.dump(overall, outfile)

    return timeline_json

def make_timeline_colours(config):

    ##make a number of hex codes that has the same number as the number of dates then assign it to config["timeline_colours"]

    config["timeline_colours"] = ['#e6959c', '#911a24']

