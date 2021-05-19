import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc



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
                sys.stderr.write(cyan(F"Error: {config["alt_seq_name"]} not found in metadata file.\n") + "Please provide a column containing alternate sequence names, or use --anonymise if you would like civet to make them for you.\n"))
                sys.exit(-1)

    elif config["anonymise"]:
        anon_dict = create_anon_ids(metadata)
        add_col_to_metadata("alternative_seq_name", anon_dict, metadata)



def create_anon_ids(metadata):

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






