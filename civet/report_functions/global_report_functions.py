
import csv
import sys
import os
import random
from civet.utils import misc
from civet.utils.log_colours import green,cyan,red

def sequence_name_parsing(report_column, anonymise, config):
    """
    parses the report group arguments 
    --report-column (Default $SEQ_NAME)
    --anonymise (Default False)
    """
    # if command line arg, overwrite config value
    misc.add_arg_to_config("report_column",report_column,config)
    misc.add_arg_to_config("anonymise",anonymise,config)

    if config["anonymise"]:
        name_dict = {}
        if "input_csv" in config:
            name_dict = create_anon_ids_from_csv(config, config["input_csv"])
        elif "ids" in config:
            name_dict = create_anon_ids_from_list(config["ids"])

        config["report_column"] = "anon_names"

        return name_dict

    elif config["report_column"]:
        if "input_csv" in config:
            with open(config["input_csv"]) as f:
                reader = csv.DictReader(f)
                if config["report_column"] not in reader.fieldnames:
                    sys.stderr.write(cyan(f"Error: {config['report_column']} not found in input metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n")
                    sys.exit(-1)
                for row in reader:
                    if '|' in row[config["report_column"]]:
                        sys.stderr.write(cyan(f"Error: {config['report_column']} contains '|' characters.\n") + "Please remove and try again.\n")
                        sys.exit(-1)
        else:
            config["report_column"] = config["background_column"]

    else:
        config["report_column"] = config["input_column"]

    return


def create_anon_ids_from_csv(config, metadata): #the names need to be swapped out

    anon_dict = {}
    name_list = []

    with open(metadata) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name_list.append(row[config["input_column"]])

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

def write_anon_names_to_file(config, name_dict):

    new_metadata = os.path.join(config["tempdir"], "query_metadata.anonymised.master.csv")

    misc.add_col_to_metadata(config["report_column"],name_dict, config["query_metadata"], new_metadata, config["input_column"], config)

    config["query_metadata"] = new_metadata

def qc_date_col(column_arg, config, metadata, metadata_name, cl_arg):

    with open(metadata) as f:
        reader = csv.DictReader(f)
    
        if config[column_arg]:
            if config[column_arg] not in reader.fieldnames:
                sys.stderr.write(cyan(f"Error: {config[column_arg]} column not found in {metadata_name} metadata file. Please specifiy which column to match with {cl_arg}`\n"))
                sys.exit(-1)

    if config[column_arg]:
        count = 0
        with open(metadata) as f:
            reader = csv.DictReader(f)
            for row in reader:
                count += 1
                if row[config[column_arg]] != "":
                    misc.check_date_format(row[config[column_arg]], count, config[column_arg])    

def parse_date_args(date_column, background_date_column, config):

    """
    parses the report group arguments:
    --date-column (default: sample_date) 
    --background-date-column (default: False, and then the same as date_column if none provided)
    """

    misc.add_arg_to_config("date_column", date_column, config)
    misc.add_arg_to_config("background_date_column", background_date_column, config)

    if "input_csv" in config:
        if not config["date_column"]:
            with open(config["input_csv"]) as f:
                reader = csv.DictReader(f)
                if "sample_date" in reader.fieldnames:
                    config["date_column"] = "sample_date"
        
    if "input_csv" in config and config["date_column"]: #so this will also scoop in "sample_date" assigned above
        qc_date_col("date_column", config, config["input_csv"], "input", "-date/--date-column")

    if not config["background_date_column"]:
        with open(config["background_csv"]) as f:
                reader = csv.DictReader(f)
                if "sample_date" in reader.fieldnames:
                    config["background_date_column"] = "sample_date"
                elif config["date_column"] and config["date_column"] in reader.fieldnames:
                    config["background_date_column"] = config["date_column"]

    if config["background_date_column"]:
        qc_date_col("background_date_column", config, config["background_csv"], "background", "-bdate/--background-date-column")


def parse_location(location_column, config):
    """
    parses the report group arguments:
    --location_column (default: country)
    """

    misc.add_arg_to_config("location_column", location_column, config)

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
     
        if config["location_column"]:
            if config["location_column"] not in headers:
                sys.stderr.write(cyan(f"Error: {config['location_column']} column not found in background metadata file. Please specify which column to match with -loc/--location`\n"))
                sys.exit(-1)
        else:
            if "country" in headers:
                config["location_column"] = "country"
