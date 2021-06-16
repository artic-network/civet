
import csv
import sys
import os
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

def qc_date_col(column_arg, config, metadata, metadata_name, cl_arg):

    with open(metadata) as f:
        reader = csv.DictReader(f)
    
        if config[column_arg]:
            if config[column_arg] not in reader.fieldnames:
                sys.stderr.write(cyan(f"Error: {column_arg} column not found in {metadata_name} metadata file. Please specifiy which column to match with {cl_arg}`\n"))
                sys.exit(-1)
        else:
            if "sample_date" in reader.fieldnames:
                config[column_arg] = "sample_date"

    if config[column_arg]:
        count = 0
        with open(metadata) as f:
            data = csv.DictReader(f)
            for line in data:
                count += 1
                if line[config[column_arg]] != "":
                    misc.check_date_format(line[config[column_arg]], count, config[column_arg])    

def parse_date_args(date_column, background_date_column, config):

    """
    parses the report group arguments:
    --date-column (default: sample_date if present, False if not)
    --background-date-column (default: sample_date if present, False if not)
    """

    misc.add_arg_to_config("date_column", date_column, config)
    misc.add_arg_to_config("background_date_column", background_date_column, config)

    if "input_csv" in config:
        qc_date_col("date_column", config, config["input_csv"], "input", "-date/--date-column")
    elif config["date_column"]:
        sys.stderr.write(cyan("You have provided a column header to find date in the input csv using '-date/--date-column' but no input csv. Please use '-bdate/--background-date-column to provide date information from the background metadata, or provide an input csv with this data\n"))
        sys.exit(-1)

    qc_date_col("background_date_column", config, config["background_csv"], "background", "-bdate/--background-date-column")


