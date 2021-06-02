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
                    new_dict["sequence_name"] = l["sequence_name"]
                    new_dict["date"] = l[col]
                    new_dict["date_type"] = col
                    dict_list.append(new_dict)
                
                
                overall[l['catchment']] = dict_list

    with open(os.path.join(config["tempdir"],'timeline_data.json'), 'w') as outfile:
        json.dump(overall, outfile)







