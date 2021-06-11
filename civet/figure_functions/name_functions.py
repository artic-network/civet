def sequence_name_qc(alt_seq_name, config):
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
