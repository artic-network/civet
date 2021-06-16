
import csv
import sys
from civet.utils import misc
from civet.utils.log_colours import green,cyan,red


def parse_and_qc_table_cols(table_content, config):

    """ 
    parses the report group argument:
    --table-content (default --background_column,--background_date_column,source,lineage,country,catchment)
    """

    processed_columns = ["source","qc_status","display_name","catchment"]
    # processed_columns = ["source","qc_status","display_name","catchment","SNP_distance","closest","SNP_list"] #update when the analysis pipeline is done

    # if command line arg, overwrite config value
    misc.add_arg_to_config("table_content",table_content,config)

    if "input_csv" in config:
        with open(config["input_csv"]) as f:
            reader = csv.DictReader(f)
            input_fieldnames = reader.fieldnames
    else:
        input_fieldnames = []

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames

    if config["table_content"]:
        if type(config["table_content"]) == str:
            content_list = config["table_content"].split(",")
        else:
            content_list = config["table_content"]
        
        for col in content_list:
            if col not in processed_columns and col not in background_fieldnames and col not in input_fieldnames:
                sys.stderr.write(cyan(f"Error: {col} column not found in metadata file\n"))
                sys.exit(-1)
           
        if config["input_column"] not in content_list: 
            content_list.insert(0,config["input_column"])

        config["table_content"] = content_list 
    else:
        sort_default_headers(input_fieldnames, background_fieldnames, config)


def sort_default_headers(input_fieldnames, background_fieldnames, config):

    basic_default_list = [config["input_column"], "lineage", "source", "catchment"]
    full_default_list = [config["background_date_column"], config["date_column"], "country"]

    header_list = basic_default_list

    for col in full_default_list:
        if col: #deals with the date columns
            if col in (background_fieldnames or col in input_fieldnames) and col not in header_list: #so if eg the date columns are the same, they don't get added twice
                header_list.append(col)

    config["table_content"] = header_list

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
 