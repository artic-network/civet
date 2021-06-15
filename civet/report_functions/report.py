import json
import csv
from mako.template import Template
from mako.lookup import TemplateLookup
import datetime as dt
from civet import __version__
from civet.report_functions import timeline_functions



def check_which_tables_produced(metadata, config): #call this in the render report

    queries_provided = False
    queries_found = False

    with open(metadata) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["source"] == "background_data":
                queries_found = True
            if row["source"] == "input_fasta":
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
            table_row = {}
            for col in config["table_content"]:
                table_row[col] = row[col]
            table_row["source"] = row["source"]
                
            query_summary_data.append(table_row)
    return query_summary_data

def make_catchment_summary_data(metadata):

    return catchment_summary_data


def get_timeline(config):

    timeline_data = report_functions.make_timeline_json(config)

    # with open(timeline_json,'r') as f:
    #     timeline_data = json.load(f)
    
    return timeline_data

def make_timeline_colours(config):

    ##make a number of hex codes that has the same number as the number of dates then assign it to config["timeline_colours"]

    return config

def make_report(metadata,report_to_generate,config):
    #all of the if statements
    #need to call this multiple times if there are multiple reports wanted

    query_summary_data = make_query_summary_data(metadata, config)
    
    catchments = [f"catchment_{i}" for i in range(1,config["catchment_count"]+1)]
    print(__version__,catchments)

    date = dt.datetime.today()

    chunks = TemplateLookup(directories=["../civet/data/report_chunks/"]) #absolute or relative works

    mytemplate = Template(filename=config["report_template"], strict_undefined=True, lookup=chunks)
    
    with open(report_to_generate, 'w') as f:
        f.write(mytemplate.render(date=date,
                                query_summary_data=query_summary_data,config=config,
                                timeline_data="",
                                catchments=catchments,
                                version=__version__)
                                )

