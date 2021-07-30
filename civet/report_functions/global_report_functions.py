
import csv
import sys
import os
import random
from civet.utils import misc
from civet.utils.log_colours import green,cyan,red
from civet.utils.config import *

def sequence_name_parsing(report_column, anonymise, config):
    """
    parses the report group arguments 
    --input-display-column (Default $SEQ_NAME)
    --anonymise (Default False)
    """
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_INPUT_DISPLAY_COLUMN,report_column,config)
    misc.add_arg_to_config(KEY_ANONYMISE,anonymise,config)

    if config[KEY_ANONYMISE]:
        name_dict = {}
        if KEY_INPUT_METADATA in config:
            name_dict = create_anon_ids_from_csv(config, config[KEY_INPUT_METADATA])
        elif KEY_IDS in config:
            name_dict = create_anon_ids_from_list(config[KEY_IDS])

        config[KEY_INPUT_DISPLAY_COLUMN] = "anon_names"

        return name_dict

    elif config[KEY_INPUT_DISPLAY_COLUMN]:
        if KEY_INPUT_METADATA in config:
            with open(config[KEY_INPUT_METADATA]) as f:
                reader = csv.DictReader(f)
                if config[KEY_INPUT_DISPLAY_COLUMN] not in reader.fieldnames:
                    sys.stderr.write(cyan(f"Error: {config[KEY_INPUT_DISPLAY_COLUMN]} not found in input metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n")
                    sys.exit(-1)
                for row in reader:
                    if '|' in row[config[KEY_INPUT_DISPLAY_COLUMN]]:
                        sys.stderr.write(cyan(f"Error: {config[KEY_INPUT_DISPLAY_COLUMN]} contains '|' characters.\n") + "Please remove and try again.\n")
                        sys.exit(-1)
        else:
            config[KEY_INPUT_DISPLAY_COLUMN] = config[KEY_BACKGROUND_ID_COLUMN]

    else:
        config[KEY_INPUT_DISPLAY_COLUMN] = config[KEY_INPUT_ID_COLUMN]

    test_duplicates = []
    if KEY_INPUT_METADATA in config:
        with open(config[KEY_INPUT_METADATA]) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row[config[KEY_INPUT_DISPLAY_COLUMN]] not in test_duplicates:
                    test_duplicates.append(row[config[KEY_INPUT_DISPLAY_COLUMN]])
                else:
                    sys.stderr.write(cyan(f"Error: {row[config[KEY_INPUT_DISPLAY_COLUMN]]} in input metadata twice. Please ensure unique IDs are used. To present these sequences together on the timeline, use the -tgc/--timeline-group-column argument\n"))
                    sys.exit(-1)
    return


def create_anon_ids_from_csv(config, metadata): #the names need to be swapped out

    anon_dict = {}
    name_list = []

    with open(metadata) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name_list.append(row[config[KEY_INPUT_ID_COLUMN]])

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

    new_metadata = os.path.join(config[KEY_TEMPDIR], "query_metadata.anonymised.master.csv")

    misc.add_col_to_metadata(config[KEY_INPUT_DISPLAY_COLUMN],name_dict, config[KEY_QUERY_METADATA], new_metadata, config[KEY_INPUT_ID_COLUMN], config)

    config[KEY_QUERY_METADATA] = new_metadata

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
    --input-date-column (default: sample_date) 
    --background-date-column (default: False, and then the same as date_column if none provided)
    """

    misc.add_arg_to_config(KEY_INPUT_DATE_COLUMN, date_column, config)
    misc.add_arg_to_config(KEY_BACKGROUND_DATE_COLUMN, background_date_column, config)

    if KEY_INPUT_METADATA in config:
        if not config[KEY_INPUT_DATE_COLUMN]:
            with open(config[KEY_INPUT_METADATA]) as f:
                reader = csv.DictReader(f)
                if "sample_date" in reader.fieldnames:
                    config[KEY_INPUT_DATE_COLUMN] = "sample_date"
        
    if KEY_INPUT_METADATA in config and config[KEY_INPUT_DATE_COLUMN]: #so this will also scoop in "sample_date" assigned above
        qc_date_col(KEY_INPUT_DATE_COLUMN, config, config[KEY_INPUT_METADATA], "input", "-idate/--input-date-column")

    if not config[KEY_BACKGROUND_DATE_COLUMN]:
        with open(config[KEY_BACKGROUND_METADATA]) as f:
                reader = csv.DictReader(f)
                if "sample_date" in reader.fieldnames:
                    config[KEY_BACKGROUND_DATE_COLUMN] = "sample_date"
                elif config[KEY_INPUT_DATE_COLUMN] and config[KEY_INPUT_DATE_COLUMN] in reader.fieldnames:
                    config[KEY_BACKGROUND_DATE_COLUMN] = config[KEY_INPUT_DATE_COLUMN]

    if config[KEY_BACKGROUND_DATE_COLUMN]:
        qc_date_col(KEY_BACKGROUND_DATE_COLUMN, config, config[KEY_BACKGROUND_METADATA], "background", "-bdate/--background-date-column")


def parse_location(location_column, config):
    """
    parses the report group arguments:
    --location_column (default: country)
    """

    misc.add_arg_to_config(KEY_BACKGROUND_LOCATION_COLUMN, location_column, config)

    with open(config[KEY_BACKGROUND_METADATA]) as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
     
        if config[KEY_BACKGROUND_LOCATION_COLUMN]:
            if config[KEY_BACKGROUND_LOCATION_COLUMN] not in headers:
                sys.stderr.write(cyan(f"Error: {config[KEY_BACKGROUND_LOCATION_COLUMN]} column not found in background metadata file. Please specify which column to match with -loc/--location`\n"))
                sys.exit(-1)
        else:
            if "country" in headers:
                config[KEY_BACKGROUND_LOCATION_COLUMN] = "country"
