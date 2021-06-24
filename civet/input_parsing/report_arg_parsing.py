import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt

from civet.report_functions import table_functions
from civet.report_functions import timeline_functions
from civet.report_functions import map_functions
from civet.report_functions import global_report_functions

def qc_report_content(config): #doesn't work with default
    reports = config["report_content"]
    to_generate = []
    to_run = []
    for report_options in reports:
        report_options = report_options.split(",")
        try:
            report_options = [int(i) for i in report_options]
        except:
            sys.stderr.write(cyan(f'Error: -rc/ --report-content should be a comma separated numerical string.\n'))
            sys.exit(-1)
        
        for i in report_options:
            if i in range(1,8):
                to_run.append(i)
            else:
                sys.stderr.write(cyan(f'Error: {i} is an invalid -rc/ --report-content option. Options 1..7 inclusive.\n'))
                sys.exit(-1)

        report_options = sorted(list(set(report_options)))
        to_generate.append(report_options)
        
    to_run = sorted(list(set(to_run)))    

    if len(to_generate)>1:
        print(green("Reports to generate:"))
        c = 0
        for i in to_generate:
            c+=1
            content_str = ",".join([str(x) for x in i])
            print(f"{c}. {content_str}")
    else:
        print(green("Report to generate:"))
        content_str = ",".join([str(x) for x in to_generate[0]])
        print(f"{content_str}")
        
    config["report_content"] = to_run
    config["reports"] = to_generate

#then at some point we need to update the treefile with these display names using jclusterfunk

def parse_global_report_options(report_content,report_column, anonymise,date_column, background_date_column,location_column, config):
    """
    parses the report group arguments 
    --report-content (Default 1,2,3)
    --report-column (Default $SEQ_NAME)
    --anonymise (Default False)
    --date-column (default: sample_date if present, False if not)
    --background-date-column (default: sample_date if present, False if not)
    --timeline-dates
    """

    # if command line arg, overwrite config value
    misc.add_arg_to_config("report_content",report_content,config)
    qc_report_content(config)

    #global report options
    name_output = global_report_functions.sequence_name_parsing(report_column, anonymise, config)
    global_report_functions.parse_date_args(date_column, background_date_column, config)
    global_report_functions.parse_location(location_column, config)

    return name_output


def parse_optional_report_content(table_content, timeline_dates, timeline_colours, config):
    #parse optional parts of report

    if 1 in config['report_content']:
        table_functions.parse_and_qc_table_cols(table_content, config)
        printable_cols = ",".join(config["table_content"])
        print(green("Metadata table will contain the following columns: ") + f"{printable_cols}\n")

    if 5 in config['report_content']:
        timeline_functions.timeline_checking(timeline_dates, config)
        timeline_functions.make_timeline_colours(timeline_colours,config)


def parse_map_options(background_map_date_restriction, background_map_column, background_map_file, query_map_file, longitude_column, latitude_column, found_in_background_data, config):

    if 6 in config['report_content']:
        map_functions.parse_background_map_options(background_map_file, background_map_date_restriction, background_map_column, found_in_background_data, config)

    if 7 in config['report_content']:
        map_functions.parse_query_map(query_map_file,longitude_column, latitude_column, found_in_background_data, config)