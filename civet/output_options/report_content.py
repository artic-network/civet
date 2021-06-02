import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt



def sequence_name_parsing(metadata, alt_seq_name, anonymise, config):
    """
    parses the report group arguments 
    --alt-seq-name (Default $SEQ_NAME)
    --anonymise (Default False)
    """
    # if command line arg, overwrite config value
    misc.add_arg_to_config("alt_seq_name",alt_seq_name,config)
    misc.add_arg_to_config("anonymise",anonymise,config)

    ##need to do some bits here - should add a column saying display name or something regardless, and then other figures can pull from that
    if config["alt_seq_name"]:
        with open(metadata) as f:
            data = csv.DictReader(f)
            if config["alt_seq_name"] not in data.fieldnames:
                sys.stderr.write(cyan(f"Error: {config['alt_seq_name']} not found in metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n")
                sys.exit(-1)

    elif config["anonymise"]:
        anon_dict = create_anon_ids(config, metadata)
        add_col_to_metadata("display_name", anon_dict, metadata) #needs updating

    else:
        add_col_to_metadata("display_name", ) #needs updating

    return config
#then at some point we need to update the treefile with these display names using jclusterfunk


def create_anon_ids(config, metadata): #the names need to be swapped out

    anon_dict = {}
    name_list = []

    with open(metadata) as f:
        data = csv.DictReader(f)
        for line in data:
            if line["query_boolean"] == "TRUE": 
                name_list.append(line[config["data_column"]])

    count = 0
    for query in random.shuffle(name_list):
        anon_dict[query] = f"sequence_{count}"

    return anon_dict


def timeline_checking(metadata, timeline_dates, config):  

    misc.add_arg_to_config("timeline_dates",timeline_dates,config) #default is None

    if config["timeline_dates"]:

        date_header_list = config["timeline_dates"].split(",")

        with open(metadata) as f:
            data = csv.DictReader(f)
            for header in date_header_list:
                if header not in data.fieldnames:
                    sys.stderr.write(cyan(f"Error: {header} (specified for use in timeline plot) not found in metadata file.\n"))
                    sys.exit(-1)

            line_count = 0
            for line in data:
                line_count += 1
                for header in date_header_list:
                    if line[header] != "":
                        misc.check_date_format(line[header], line_count, header)


    return config


def make_timeline_json(config):

    date_cols = config["timeline_dates"].split(",")

    overall = defaultdict(dict)

    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['query_boolean'] == "TRUE":
                            
                if l['catchment'] in overall:
                    dict_list = overall[l['catchment']]
                else:
                    dict_list = []
                
                for col in date_cols:
                    new_dict = {}
                    new_dict["sequence_name"] = l["display_name"]
                    new_dict["date"] = l[col]
                    new_dict["date_type"] = col
                    dict_list.append(new_dict)
                
                
                overall[l['catchment']] = dict_list

    with open(os.path.join(config["tempdir"],'timeline_data.json'), 'w') as outfile:
        json.dump(overall, outfile)

def parse_date_columns(metadata, data_date_column, input_date_column, config): #needs to be done on combined metadata file currently - can be edited to be done earlier on if the combined metadata is done too late

    misc.add_arg_to_config("data_date_column", data_date_column, config)
    misc.add_arg_to_config("input_date_column", input_date_column, config)

    with open(metadata) as f:
        reader = csv.DictReader(f)
        
        if data_date_column not in reader.fieldnames:
            sys.stderr.write(cyan(f"Error: {data_date_column} column not found in metadata file. Please specifiy which column to match with `-ddcol/--data-date-column.`\n"))
            sys.exit(-1)
        elif input_date_column not in reader.fieldnames:
            sys.stderr.write(cyan(f"Error: {input_date_column} column not found in metadata sfile. Please specifiy which column to match with `-idatcol/--input-date-column.`\n"))
            sys.exit(-1)



def parse_and_qc_table_cols(metadata, found_table_cols, provided_table_cols, config): #needs to be done on combined metadata file

    """
    parses the report group arguments:
    --found-table-cols (Default --data_column,--data_date_column,lineage,country,catchment)
    --provided-table-cols (Default --fasta_column,--input_date_column,closest,SNP_distance,SNP_list")
    """

    # if command line arg, overwrite config value
    misc.add_arg_to_config("found_table_cols",found_table_cols,config)
    misc.add_arg_to_config("provided_table_cols",provided_table_cols,config)

    with open(metadata) as f:
        reader = csv.DictReader()
        fieldnames = reader.fieldnames
    
    if config["found_table_cols"]:
        for col in config["found_table_cols"]:
            if col != "catchment":
                if col not in fieldnames:
                    sys.stderr.write(cyan(f"Error: {col} column not found in metadata file. Please specifiy which column to match with `-ddcol/--data-date-column.`\n"))
                    sys.exit(-1)
    else:
        config["found_table_cols"] = ",".join([config["data_column"], config["data_date_column"],"lineage","country","catchment"])
    
    #don't need to check this - the others are already checked, and then the other stuff is within the processed metadata 
    #think about this - maybe don't allow them to do this? not sure
    if not config["provided_table_cols"]: 
        config["provided_table_cols"] = ",".join([config["fasta_column"], config['input_date_column'],"closest","SNP_distance","SNP_list"])








