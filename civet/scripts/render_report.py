#things that need to be specified:

#plausibly can pass the config to the mako template

#data_column
#date of report ($date)
#list of rows for each table (query_table_cols, background_table_cols) DONE






def check_which_tables_produced(metadata, config): #call this in the render report

    queries_provided = False
    queries_found = False

    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            if l["source"] == "background_data":
                queries_found = True
            if l["source"] == "input_fasta":
                queries_provided = True

    config["queries_found"] = queries_found
    config["queries_provided"] = queries_provided

def process_catchments():

    #output is list of catchment IDs

    #dict of 
    #catchment_id: catchment_1
    #treeString: newick_string


