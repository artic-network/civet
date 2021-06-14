
import csv
import sys
from civet.utils import misc
from civet.utils.log_colours import green,cyan,red


def parse_and_qc_table_cols(found_seq_table, novel_seq_table, config): 
    """
    parses the report group arguments:
    --found-seq-table (Default --data_column,--background_date_column,lineage,country,catchment)
    --novel-seq-table (Default --fasta_column,--date_column,closest,SNP_distance,SNP_list")
    """

    processed_columns = ["source","qc_status","display_name","catchment"]
    # processed_columns = ["source","qc_status","display_name","catchment","SNP_distance","closest","SNP_list"] #update when the analysis pipeline is done

    # if command line arg, overwrite config value
    misc.add_arg_to_config("found_seq_table",found_seq_table,config)
    misc.add_arg_to_config("novel_seq_table",novel_seq_table,config)

    if "input_csv" in config:
        with open(config["input_csv"]) as f:
            reader = csv.DictReader(f)
            input_fieldnames = reader.fieldnames
    else:
        input_fieldnames = []

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames
    
    if config["found_seq_table"]:
        for col in config["found_seq_table"]:
            if col not in processed_columns and col not in background_fieldnames and col not in input_fieldnames:
                sys.stderr.write(cyan(f"Error: {col} column not found in metadata file\n"))
                sys.exit(-1)
    else:
        if config["background_date_column"]:
            config["found_seq_table"] = ",".join([config["data_column"], config["background_date_column"],"catchment"])
        else:
            config["found_seq_table"] = ",".join([config["data_column"],"catchment"])
    

    if config["novel_seq_table"]:
        for col in config["novel_seq_table"]:
            if col not in processed_columns and col not in input_fieldnames:
                sys.stderr.write(cyan(f"Error: {col} column not found in metadata file\n"))
                sys.exit(-1)

    else: 
        # config["novel_seq_table"] = ",".join([config["fasta_column"], config['date_column'],"closest","catchment","SNP_distance","SNP_list"])
        if config["date_column"]:
            config["novel_seq_table"] = ",".join([config["fasta_column"], config['date_column'],"catchment"])
        else:
            config["novel_seq_table"] = ",".join([config["fasta_column"],"catchment"])

def parse_input_date(date_column, config):

    misc.add_arg_to_config("date_column", date_column, config)
    if "input_csv" in config:

        with open(config["input_csv"]) as f:
            reader = csv.DictReader(f)
        
            if config["date_column"]:
                if config["date_column"] not in reader.fieldnames:
                    sys.stderr.write(cyan(f"Error: {date_column} column not found in query metadata file. Please specifiy which column to match with `-date/--date-column.`\n"))
                    sys.exit(-1)
            else:
                if "sample_date" in reader.fieldnames:
                    config["date_column"] = "sample_date"

        if config["date_column"]:
            count = 0
            with open(config["input_csv"]) as f:
                data = csv.DictReader(f)
                for line in data:
                    count += 1
                    if line[config["date_column"]] != "":
                        misc.check_date_format(line[config["date_column"]], count, config["date_column"])    
    elif config["date_column"]:
        sys.stderr.write(cyan("You have provided a column header to find date in the input csv using '-date/--date-column' but no input csv. Please use '-bdate/--background-date-column to provide date information from the background metadata.\n"))
        sys.exit(-1)

def parse_background_date(background_date_column, config):

    misc.add_arg_to_config("background_date_column", background_date_column, config)

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
    
        if config["background_date_column"]:
            if config["background_date_column"] not in reader.fieldnames:
                sys.stderr.write(cyan(f"Error: {background_date_column} column not found in background metadata file. Please specifiy which column to match with `-bdate/--background-date-column.`\n"))
                sys.exit(-1)
        else:
            if "sample_date" in reader.fieldnames:
                config["background_date_column"] = "sample_date"

    if config["background_date_column"]:
        count = 0
        with open(config["background_csv"]) as f:
            data = csv.DictReader(f)
            for line in data:
                count += 1
                if line[config["background_date_column"]] != "":
                    misc.check_date_format(line[config["background_date_column"]], count, config["background_date_column"])    