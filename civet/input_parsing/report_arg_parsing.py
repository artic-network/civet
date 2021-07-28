import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt

from civet.report_functions import table_functions
from civet.report_functions import timeline_functions
from civet.report_functions import map_functions
from civet.report_functions import global_report_functions
from civet.utils.config import *

def parse_preset_options(config):
    report_config = []

    valid_preset = {
        "the_usual":"1,2,3,4,5", 
        "the_works":"1,2,3,4,5,7", 
        "the_whole_shebang":"1,2,3,4,5,6,7", 
        "hold_the_sauce":"1,2"
        }
    report_preset = config[KEY_REPORT_PRESET]
    if not type(report_preset) == list:
        report_preset = report_preset.split(" ")
    
    for preset_option in report_preset:
        if preset_option not in valid_preset:
            preset_str = "- "
            for i in valid_preset:
                preset_str+=f"{i} ({','.join(valid_preset[i])})\n"
            sys.stderr.write(cyan(f'Error: valid `-rp/--report-preset` options are\n')+
            """
- the_usual (1,2,3,4,5)
- the_works (1,2,3,4,5,7) 
- the_whole_shebang (1,2,3,4,5,6,7)
- hold_the_sauce (1,2)  
            """)
            sys.exit(-1)
        else:
            report_config.append(valid_preset[preset_option])
    config[KEY_REPORT_CONTENT] = report_config

def qc_report_content(config):
    reports = config[KEY_REPORT_CONTENT]
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
    pretty_options = find_pretty_report_options()

    if len(to_generate)>1:
        print(green("Reports to generate:"))
        c = 0
        for i in to_generate:
            c+=1
            content_str = ", ".join([pretty_options[str(x)] for x in i])
            print(f"{c}. {content_str}")
    else:
        print(green("Report to generate:"))
        content_str = ", ".join([pretty_options[str(x)] for x in to_generate[0]])
        print(f"{content_str}")
        
    config[KEY_REPORT_CONTENT] = to_run
    config[KEY_REPORTS] = to_generate

#then at some point we need to update the treefile with these display names using jclusterfunk

def parse_global_report_options(report_content,report_preset,report_column, anonymise,date_column, background_date_column,location_column, config):
    """
    parses the report group arguments 
    --report-content (Default 1,2,3)
    --input-display-column (Default $SEQ_NAME)
    --anonymise (Default False)
    --input-date-column (default: sample_date if present, False if not)
    --background-date-column (default: sample_date if present, False if not)
    --timeline-dates
    """

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_REPORT_CONTENT,report_content,config)
    misc.add_arg_to_config(KEY_REPORT_PRESET,report_preset,config)
    if config[KEY_REPORT_PRESET]:
        parse_preset_options(config)
    qc_report_content(config)

    #global report options
    name_output = global_report_functions.sequence_name_parsing(report_column, anonymise, config)
    global_report_functions.parse_date_args(date_column, background_date_column, config)
    global_report_functions.parse_location(location_column, config)

    return name_output

def is_hex(code):
    hex_code = False
    if code.startswith("#") and len(code)==7:
        values = code[1:].upper()
        for value in values:
            if value in "ABCDEF0123456789":
                hex_code = True
            else:
                hex_code = False
                break
                
    return hex_code

def get_acceptable_colours(config):

    acceptable_colours = []

    with open(config["html_colours"]) as f:
        reader = csv.DictReader(f)
        for row in reader:
            acceptable_colours.append(row["name"].lower())
            # acceptable_colours.append(row["hex"].lower())

    return acceptable_colours

def colour_checking(config):
    acceptable_colours = get_acceptable_colours(config)
    cmap = config["colour_map"]
    if not type(cmap) == "list":
        if ',' in cmap:
            cmap = cmap.split(",")
            for colour in cmap:
                if not is_hex(colour) and colour.lower() not in acceptable_colours:
                    sys.stderr(cyan(f"Invalid colour code: ") + f"{colour}\n")
                    sys.exit(-1)
        else:
            cmap = cmap.split(",")
            if not is_hex(cmap[0]) and colour.lower() not in acceptable_colours:
                sys.stderr(cyan(f"Invalid colour code: ") + f"{cmap[0]}\nPlease provide a comma-separated string of HEX codes or see htmlcolorcodes.com/color-names for `-cmap/--colour-map`.\n")
                sys.exit(-1)
    config["colour_map"] = cmap
    if not is_hex(config["colour_theme"]) and config["colour_theme"].lower() not in acceptable_colours:
        sys.stderr(cyan(f"Invalid HEX colour code: ") + f"{config['colour_theme']}\nPlease provide a valid HEX code or see htmlcolorcodes.com/color-names for `-ct/--colour-theme`.\n")
        sys.exit(-1)

def parse_tree_options(tree_annotations,max_tree_size, config):
    misc.add_arg_to_config(KEY_TREE_ANNOTATIONS,tree_annotations, config)
    misc.add_arg_to_config("max_tree_size",max_tree_size, config)

    try:
        config["max_tree_size"] = int(config["max_tree_size"])
    except:
        sys.stderr(cyan(f"Error: `-mq/--max-tree-size` must be an integer."))
        sys.exit(-1)

    if not type(config[KEY_TREE_ANNOTATIONS])==list:
        config[KEY_TREE_ANNOTATIONS] = config[KEY_TREE_ANNOTATIONS].split(',')

    for col in config[KEY_TREE_ANNOTATIONS]:
        if col not in config["query_csv_header"] and col not in config[KEY_MUTATIONS]:
            sys.stderr(cyan(f"Error: `{col}`` column not provided for tree annotations."))
            sys.exit(-1)
    
    config[KEY_TREE_ANNOTATIONS] = " ".join(config[KEY_TREE_ANNOTATIONS])


def parse_optional_report_content(table_content,mutations, timeline_dates, timeline_group_column, colour_theme, colour_map, config):
    #parse optional parts of report
    

    misc.add_arg_to_config("colour_theme",colour_theme,config)
    misc.add_arg_to_config("colour_map",colour_map,config)

    colour_checking(config)

    if 1 in config[KEY_REPORT_CONTENT]:
        table_functions.parse_and_qc_table_cols(table_content,mutations, config)
        printable_cols = ",".join(config[KEY_QUERY_TABLE_CONTENT])
        print(green("Metadata table will contain the following columns: ") + f"{printable_cols}")

    if 5 in config[KEY_REPORT_CONTENT]:
        timeline_functions.timeline_checking(timeline_dates, timeline_group_column, config)


def parse_map_options(background_map_date_range, background_map_column, background_map_file, centroid_file, background_map_location, query_map_file, longitude_column, latitude_column, found_in_background_data, config):

    if 6 in config[KEY_REPORT_CONTENT]:
        map_functions.parse_background_map_options(background_map_file, centroid_file, background_map_date_range, background_map_column, background_map_location, found_in_background_data, config)

    if 7 in config[KEY_REPORT_CONTENT]:
        map_functions.parse_query_map(query_map_file,longitude_column, latitude_column, found_in_background_data, config)


def find_pretty_report_options():
    option_dict = {}
    option_dict["1"] = "Query tables"
    option_dict["2"] = "catchment tables"
    option_dict["3"] = "trees"
    option_dict["4"] = "snipit"
    option_dict["5"] = "timelines"
    option_dict["6"] = "background diversity maps"
    option_dict["7"] = "query maps"

    return option_dict