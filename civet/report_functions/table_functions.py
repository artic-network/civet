import csv
import sys
from civet.utils import misc
from civet.utils.log_colours import green,cyan,red
from civet.utils.config import *

def parse_mutations(config):
    if config[KEY_MUTATIONS]:
        if not type(config[KEY_MUTATIONS]) == list:
            config[KEY_MUTATIONS] = config[KEY_MUTATIONS].split(",")

        for mutation in config[KEY_MUTATIONS]:
            if not ':' in mutation:
                sys.stderr.write(cyan(f"Error: invalid mutation specified {mutation}\n"))
                sys.exit(-1)
            else:
                config[KEY_QUERY_TABLE_CONTENT].append(mutation)

def parse_and_qc_table_cols(table_content,mutations, config):

    """ 
    parses the report group argument:
    --query-table-content (default --background_column,--background_date_column,source,lineage,country,catchment)
    """
    
    processed_columns = ["source","qc_status",config[KEY_INPUT_DISPLAY_COLUMN],"catchment"]
    # processed_columns = ["source","qc_status","config[KEY_INPUT_DISPLAY_COLUMN]","catchment","SNP_distance","closest","SNP_list"] #update when the analysis pipeline is done

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_QUERY_TABLE_CONTENT,table_content,config)
    misc.add_arg_to_config(KEY_MUTATIONS,mutations,config)

    if KEY_INPUT_METADATA in config:
        with open(config[KEY_INPUT_METADATA]) as f:
            reader = csv.DictReader(f)
            input_fieldnames = reader.fieldnames
    else:
        input_fieldnames = []

    with open(config[KEY_BACKGROUND_METADATA]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames

    if config[KEY_QUERY_TABLE_CONTENT]:
        if type(config[KEY_QUERY_TABLE_CONTENT]) == str:
            content_list = config[KEY_QUERY_TABLE_CONTENT].split(",")
        else:
            content_list = config[KEY_QUERY_TABLE_CONTENT]
        
        for col in content_list:
            if col not in processed_columns and col not in background_fieldnames and col not in input_fieldnames:
                sys.stderr.write(cyan(f"Error: {col} column not found in metadata file\n"))
                sys.exit(-1)
           
        if config[KEY_INPUT_DISPLAY_COLUMN] not in content_list:
            content_list.insert(0,config[KEY_INPUT_DISPLAY_COLUMN])

        config[KEY_QUERY_TABLE_CONTENT] = content_list 
    else:
        sort_default_headers(input_fieldnames, background_fieldnames, config)

    parse_mutations(config)

    config[KEY_FASTA_TABLE_CONTENT] = [config[KEY_INPUT_DISPLAY_COLUMN],"seq_N_content","seq_length"]


def sort_default_headers(input_fieldnames, background_fieldnames, config):

    
    basic_default_list = [config[KEY_INPUT_DISPLAY_COLUMN], "lineage", "source", "catchment"]
    full_default_list = [config[KEY_BACKGROUND_DATE_COLUMN], config[KEY_INPUT_DATE_COLUMN], "country", "adm1", "suggested_adm2_grouping", "adm2"]

    header_list = basic_default_list

    for col in full_default_list:
        if col: #deals with the date columns
            if col in (background_fieldnames or col in input_fieldnames) and col not in header_list: #so if eg the date columns are the same, they don't get added twice
                header_list.append(col)

    config[KEY_QUERY_TABLE_CONTENT] = header_list







