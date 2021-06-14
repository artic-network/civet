import csv
import sys
import os
import random
from civet.utils import misc
from civet.utils.log_colours import green,cyan


def sequence_name_parsing(report_column, anonymise, config):
    """
    parses the report group arguments 
    --report-column (Default $SEQ_NAME)
    --anonymise (Default False)
    """
    # if command line arg, overwrite config value
    misc.add_arg_to_config("report_column",report_column,config)
    misc.add_arg_to_config("anonymise",anonymise,config)

    name_dict = {}

    if "input_csv" in config:
        metadata = config["input_csv"] 
        if config["anonymise"]:
            name_dict = create_anon_ids_from_csv(config, metadata)

        elif config["report_column"]:
            with open(metadata) as f:
                data = csv.DictReader(f)
                if config["report_column"] not in data.fieldnames:
                    sys.stderr.write(cyan(f"Error: {config['report_column']} not found in input metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n")
                    sys.exit(-1)
                else:
                    for l in data:
                        name_dict[l[config['input_column']]] = l[config['report_column']]

        else:
            with open(metadata) as f: 
                data = csv.DictReader(f)
                for l in data:
                    name_dict[l[config['input_column']]] = l[config['input_column']]


    elif "ids" in config:
        if config["anonymise"]:
            name_dict = create_anon_ids_from_list(config["ids"])
        else:
            for name in config["ids"]:
                name_dict[name] = name


    return config, name_dict


def create_anon_ids_from_csv(config, metadata): #the names need to be swapped out

    anon_dict = {}
    name_list = []

    with open(metadata) as f:
        data = csv.DictReader(f)
        for line in data:
            name_list.append(line[config["input_column"]])

    random.shuffle(name_list)

    count = 0
    for query in name_list:
        count += 1
        anon_dict[query] = f"sequence_{count}"

    return anon_dict

def create_anon_ids_from_list(input_list):

    random.shuffle(input_list)

    count = 0
    for query in input_list:
        count += 1
        anon_dict[query] = f"sequence_{count}"

    return anon_dict

def write_names_to_file(config, name_dict):

    new_metadata = os.path.join(config["tempdir"], "query_metadata.name_procssed.master.csv")

    misc.add_col_to_metadata("display_name",name_dict, config["query_metadata"], new_metadata, config["input_column"], config)

    config["query_metadata"] = new_metadata
