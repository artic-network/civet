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

    if config["alt_seq_name"]:
        with open(metadata) as f:
            data = csv.DictReader(f)
            if config["alt_seq_name"] not in data.fieldnames:
                sys.stderr.write(cyan(f"Error: {config['alt_seq_name']} not found in metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n")
                sys.exit(-1)

    elif config["anonymise"]:
        anon_dict = create_anon_ids(metadata)
        add_col_to_metadata("alternative_seq_name", anon_dict, metadata)

    return config



def create_anon_ids(metadata): #the names need to be swapped out

    anon_dict = {}
    name_list = []

    with open(metadata) as f:
        data = csv.DictReader(f)
        for line in data:
            if line["query_status"]: #not sure what the header/format is
                name_list.append(line["sequence_name"]) #or whatever

    count = 0
    for query in random.shuffle(name_list):
        anon_dict[query] = f"sequence_{count}"

    return anon_dict

def check_date_format(to_check, line_count, header):

    date_format = "%Y-%m-%d"
    
    try:
        dt.datetime.strptime(to_check, date_format).date()
    except:
        print(qcfunk.cyan(f'date {to_check} on line {line_count} in column {header} in incorrect format. Please use YYYY-MM-DD'))
        sys.exit(-1)

def timeline_checking(metadata, timeline_dates, config):  

    misc.add_arg_to_config("timeline_dates",timeline_dates,config)

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
                        check_date_format(line[header], line_count, header)


    return config










