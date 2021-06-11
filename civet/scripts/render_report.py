import json
import csv
from mako import Template
from mako.lookup import TemplateLookup
import datetime as dt



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

    return catchment_list, tree_strings

def make_query_summary_data(metadata, config):
    query_summary_data = []
    with open(metadata, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cluster_dict = {}
            for col in config["provided_table_cols"]:
                cluster_dict[col] = row[col]
            for col in config["found_table_cols"]:
                cluster_dict[col] = row[col]
            cluster_dict["source"] = row["source"]
                
            query_summary_data.append(cluster_dict)
    return query_summary_data

def make_catchment_summary_data(metadata):

    return catchment_summary_data


def get_timeline(timeline_json):
    with open(timeline_json,'r') as f:
        timeline_data = json.load(f)
    return timeline_data

def make_timeline_colours(config):

    ##make a number of hex codes that has the same number as the number of dates then assign it to config["timeline_colours"]

    return config

def make_report():
#all of the if statements
#need to call this multiple times if there are multiple reports wanted
    query_summary_data = make_summary_data(metadata, config)
    timeline_data = get_timeline(timeline_json) 
    config = make_timeline_colours(config)
    catchment_lst, tree_strings = process_catchments()

    date = dt.datetime.today()

    # mylookup = TemplateLookup(directories=["../civet/data/report_chunks/"]) #absolute or relative works

    # mytemplate = Template(filename="test_template.mako", strict_undefined=True, lookup=mylookup)
    # f = open(f"test_report.html", 'w')
    # f.write(mytemplate.render(date=date,query_summary_data=query_summary_data,config=config,
    #                           timeline_data = timeline_data,
    #                           catchments=["catchment_1","catchment_2"],
    #                           version="3.0"))
    # f.close()

