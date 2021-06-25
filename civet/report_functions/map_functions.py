import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
import json
import datetime as dt

def parse_query_map(query_map_file, longitude_column, latitude_column, found_in_background_data, config):
    """
    parses map group arguments:
    --longitude_column (default=longitude)
    --latitude_column (default=latitude)
    """

    misc.add_arg_to_config("longitude_column",longitude_column,config)
    misc.add_arg_to_config("latitude_column",latitude_column, config)

    if "input_csv" in config:
        with open(config["input_csv"]) as f:
            reader = csv.DictReader(f)
            input_fieldnames = reader.fieldnames
    else:
        input_fieldnames = []

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames

    long_input = False
    lat_input = False
    if config['longitude_column'] in input_fieldnames: #this is done for efficiency later on as bool checks are faster than looping through a list
        long_input = True 
    if config['latitude_column'] in input_fieldnames:
        lat_input = True

    if not long_input and config["longitude_column"] not in background_fieldnames:
        sys.stderr.write(cyan(f"Error: {config['longitude_column']} not found in either metadata file for mapping queries\n") + "\n")
        sys.exit(-1)
    elif not lat_input and config["latitude_column"] not in background_fieldnames:
        sys.stderr.write(cyan(f"Error: {config['latitude_column']} not found in either metadata file for mapping queries\n") + "\n")
        sys.exit(-1) 
    
    else:
        #checking that there is some data in either the input or the background metadata in those columns
        seq_to_longlat_check = {}
        if config['longitude_column'] in background_fieldnames:
            for query, metadata in found_in_background_data.items():
                if query not in seq_to_longlat_check:
                    seq_to_longlat_check[query] = [0,0]
                if metadata[config['longitude_column']] != "":
                    seq_to_longlat_check[query][0] = 1
                     
        if config['latitude_column'] in background_fieldnames:
            for query, metadata in found_in_background_data.items():
                if query not in seq_to_longlat_check:
                    seq_to_longlat_check[query] = [0,0]
                if metadata[config['latitude_column']] != "":
                    seq_to_longlat_check[query][1] = 1
        
        if long_input or lat_input:
            with open(config["input_csv"]) as f:
                data = csv.DictReader(f)
                for line in data:
                    if line[config["input_column"]] not in seq_to_longlat_check:
                        seq_to_longlat_check[line[config["input_column"]]] = [0,0]
                    if long_input:
                        if line[config['longitude_column']] != "":
                            seq_to_longlat_check[line[config["input_column"]]][0] = 1
                    if lat_input:
                        if line[config["latitude_column"]] != "":
                            seq_to_longlat_check[line[config["input_column"]]][1] = 1

        data_present = False
        for seq, data in seq_to_longlat_check.items():
            if sum(data) == 2:
                data_present = True
                break
        if not data_present:
            sys.stderr.write(cyan(f"Error: no query with longitude and latitude information contained in the columns provided ({config['longitude_column'], config['latitude_column']}), so map cannot be produced.\n"))
            sys.exit(-1)

        if config["query_map_file"]:
            parse_map_file_arg(query_map_file, "query_map_file", config) #don't need to QC this because it doesn't matter what's in it if only the query map
        else:
            if config["civet_mode"] == "CLIMB":
                map_file = "uk_map.json"
            else:
                map_file = "adm0_global.json"
            
            config["query_map_file"] = map_file


def parse_background_map_column(background_map_column, config):

    misc.add_arg_to_config("background_map_column", background_map_column, config)  

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames
        
    if config["background_map_column"]:
        if config["background_map_column"] not in background_fieldnames: 
            sys.stderr.write(cyan(f"Error: {config['background_map_column']} not found in background metadata file for mapping background diversity\n") + "\n")
            sys.exit(-1)
    
    else:
        if config["civet_mode"] == "CLIMB" and "suggested_adm2_grouping" in background_fieldnames: 
            config["background_map_column"] = "suggested_adm2_grouping"
        elif "adm1" in background_fieldnames:
            config["background_map_column"] = "adm1"
        elif config["location_column"]:
            config["background_map_column"] = config["location_column"]
        else:
            sys.stderr.write(cyan(f"Error: no field found in background metadata file for mapping background diversity. Please provide one with -mapcol/--background-map-column.\n") + "\n")
            sys.exit(-1)

def parse_date_restriction(background_map_date_restriction, config):

    misc.add_arg_to_config("background_map_date_restriction", background_map_date_restriction, config)   

    if config["background_map_date_restriction"]:
        if config["background_date_column"]:
            date_range = config["background_map_date_restriction"].split(":")
            if len(date_range) > 1:
                try:
                    start_date = dt.datetime.strptime(date_range[0], "%Y-%m-%d").date()
                except:
                    sys.stderr.write(cyan("Start date in background map date restriction in incorrect format. Please use YYYY-MM-DD"))
                    sys.exit(-1)
                try:
                    end_date = dt.datetime.strptime(date_range[1], "%Y-%m-%d").date()
                except:
                    sys.stderr.write(cyan("End date in background map date restriction in incorrect format. Please use YYYY-MM-DD"))
                    sys.exit(-1)

            else:
                if not date_range[0].isdigit():
                    sys.stderr.write(cyan("Date window must be a single integer. If you have attempted provided a date, you must provide both the start and end date separated by a colon."))
                    sys.exit(-1)
                
                else:
                    start_date, end_date = do_date_window(date_range[0], found_in_background_metadata, config)

            config["start_date"] = start_date
            config["end_date"] = end_date

            print(green(f"Restricting background map analysis between the dates {start_date} and {end_date}"))

            with open(config["background_csv"]) as f:
                data = csv.DictReader(f)
                data_count = 0
                enough_data = False
                for line in data:
                    if line[config["background_date_column"]] != "":
                        date = dt.datetime.strptime(line[config["background_date_column"]], "%Y-%m-%d").date()
                        if date >= start_date and date <= end_date and line[config["background_map_column"]] != "":
                            data_count += 1
                        if data_count > 50:
                            enough_data = True
                            break
                if not enough_data:
                    sys.stderr.write(cyan(f"Error: fewer than 50 sequences in column {config['background_map_column']} between dates {start_date} and {end_date}. If you want to map background diversity, consider using a different geographical scale using '-maploc/--background-map-location' or a larger time period.\n"))
                    sys.exit(-1)
        else:
            sys.stderr.write(cyan(f"Error: Date restriction defined for background lineage diversity mapping, but no date column provided. Please use '-bdate/--background-date-column' to specify this column. \n") + "\n")
            sys.exit(-1)
    else:
        print(green("No date restriction provided, analysing background lineage diversity for the whole pandemic.\n"))
        config["start_date"] = dt.datetime(2019,12,1).date()
        config["end_date"] = dt.datetime.today().date()

def parse_background_map_options(background_map_file, background_map_date_restriction, background_map_column, found_in_background_metadata, config):

    """
    parses map group arguments:
    --background-map-date-restriction (default=no restriction)
    --background-map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode == CLIMB)
    # --map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode==CLIMB)
    """

    parse_background_map_column(background_map_column, config)
    parse_date_restriction(background_map_date_restriction, config)    

    if config["background_map_file"]:
        parse_map_file(background_map_file,"background_map_file", config)
    
    map_file = qc_map_file_for_background_map(config)

    config["background_map_file"] = map_file


def do_date_window(date_window, found_in_background_metadata, config):

    start_date = None
    end_date = None

    date_diff = dt.timedelta(int(date_window))
        
    dates = set()
    if "input_csv" in config and config['date_column']:
        with open(config["input_csv"]) as f:
            data = csv.DictReader(f)
            for l in data:
                if l[config['date_column']] != "": 
                    dates.add(dt.datetime.strptime(l[config['date_column']], '%Y-%m-%d').date())

    for seq, metadata in found_in_background_metadata.items(): 
        if metadata[config["background_date_column"]] != "":
            dates.add(dt.datetime.strptime(metadata[config["background_date_column"]],'%Y-%m-%d').date())

    if len(dates) == 0:
        sys.stderr.write(cyan("Error: no dates in queries to restrict using date window. Please provide absolute dates in the format 'YYYY-MM-DD:YYYY-MM-DD' using the argument '-daterestric/--background-map-date-restriction'\n"))
        sys.exit(-1)
          
    earliest = sorted(dates)[0]
    latest = sorted(dates)[-1]

    start_date = earliest - date_diff
    end_date = latest + date_diff

    return start_date, end_date

def parse_map_file_arg(map_file_arg, arg_name, config): 

    """
    parses map group arguments:
    --background-map-file (default=UK_map if civet_mode == CLIMB, relevant level of world_map (ie adm0 or adm1) if not)
    --query-map-file
    """
    misc.add_arg_to_config(arg_name, map_file_arg, config)

    if config[arg_name]:
        if not os.path.exists(config[arg_name]):
            sys.stderr.write(cyan(f"{config[arg_name]} cannot be found. Please check the path and try again.\n"))
            sys.exit(-1)
        if not config[arg_name].endswith(".json") and not config[arg_name].endswith(".geojson"):
            sys.stderr.write(csyan(f"{config[arg_name]} must be in the format of a geojson. You can use mapshaper.org to convert between file formats.\n"))
            sys.exit(-1)
 

    return map_file

def qc_map_file_for_background_map(config):

    if config["verbose"]:
        print("Beginning checks for background map")

    #winding in the geojsons

    if config["background_map_file"]:

        with open(config["background_map_file"]) as f:
            geodata = json.load(f)
        
        headers = geodata["features"][0]["properties"].keys()
        
        if config["background_map_column"] not in headers:
            sys.stderr.write(cyan(f"{config['background_map_column']} not found in custom shapefile.\n"))
            sys.exit(-1)
        else:
            acceptable_locations = []
            for item in geodata['features']:
                acceptable_locations.append(item["properties"][config["background_map_column"]])
    
    else:
        if config["civet_mode"] == "CLIMB":
            map_file = "uk_map.json"
            uk_cols = ["suggested_adm2_grouping", "adm1", "adm2"]
            if config["background_map_column"] not in uk_cols:
                sys.stderr.write(cyan(f'{config["background_map_column"]} not in default UK map file.\n Options allowed are "suggested_adm2_grouping","adm1" or "adm2".  Alternatively, please provide a custom geojson containing this column to use it using --map-file\n'))
                sys.exit(-1)
        else: 
            if config["background_map_column"] == "adm1":
                map_file = "adm1_global.json"
            elif config["background_map_column"] == "country" or config["background_map_column"] == "adm0" or config["background_map_column"] == "ISO":
                map_file = "adm0_global.json"
            else:
                sys.stderr.write(cyan(f"{config['background_map_column']} not in default map file. Please use country/adm0 or adm1 or provide your own shape file using --map-file\n"))
                sys.exit(-1)

        acceptable_locations = get_acceptable_locations(map_file, config)

    check_set = set()
    with open(config["background_csv"]) as f:
        data = csv.DictReader(f)
        for line in data:
            location_value = line[config["background_map_column"]]
            if config["background_date_column"]:
                date_value = line[config["background_date_column"]]
                if date_value != "":
                    date = dt.datetime.strptime(date_value, "%Y-%m-%d").date()
                    if date >= config["start_date"] and date <= config["end_date"]:
                        if location_value != "":
                            if config["civet_mode"] == "CLIMB": 
                                if line["country"] == "UK" and location_value != "Needs_manual_curation" and "|" not in location_value:
                                    check_set.add(location_value)
                            else:
                                check_set.add(location_value)
            else:
                if config["civet_mode"] == "CLIMB": 
                    if line["country"] == "UK" and location_value != "Needs_manual_curation" and "|" not in location_value:
                        check_set.add(location_value)
                else:
                    check_set.add(location_value)

    for loc in check_set: #this needs to be different if it's in the UK, only check UK locations
        if loc not in acceptable_locations:
            if config["background_map_file"]:
                sys.stderr.write(f'{loc} is an invalid location. Please ensure that the metadata values match up to the map file you have provided.\n')
                sys.exit(-1)
            elif config['civet_mode'] == "CLIMB":
                sys.stderr.write(f'{loc} is an invalid location. If you are using the default background metadata, please contact Verity Hill (verity.hill@ed.ac.uk)\n')
                sys.exit(-1)
            else:
                sys.stderr.write(f"{loc} isn't in our map file. Please see a list of currently accepted locations here: [link]. If you can't find your country's data on that list, please open a github issue on the civet repo and we will get to it as soon as we can.\n")
                sys.exit(-1)

    if config["verbose"]:
        print("Finished with checks for background map")

    return map_file

def get_acceptable_locations(map_file, config):

    if config["civet_mode"] == "CLIMB":
        if config["background_map_column"] == "adm1":
            acceptable_locations = ["Scotland", "Wales", "Northern_Ireland", "England", "Jersey", "Guernsey", "Isle_of_Man", "Falkland_Islands", "Gibraltar"]
        else:
            acceptable_locations = []
            with open(config["uk_acceptable_values"]) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    acceptable_locations.append(row[config["background_map_column"]])

    else:
        acceptable_locations = set()
        with open(config["global_acceptable_values"]) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acceptable_locations.add(row[config["background_map_column"]])
    
    return acceptable_locations
                

def make_query_map_json(config, metadata): 
#colour by is dynamic
    lat_col = config["latitude_column"]
    long_col = config["longitude_column"]
    name_col = config["background_column"] 

    overall_dict = defaultdict(dict)

    all_queries = {}
    all_queries['queries'] = []

    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            per_seq_dict = {}
            per_seq_dict["sequence_name"] = l[name_col]
            per_seq_dict["latitude"] = l[lat_col]
            per_seq_dict["longitude"] = l[long_col]
            
            all_queries['queries'].append(per_seq_dict)

    with open(os.path.join(config["tempdir"],'query_map_data.json'), 'w') as outfile:
        json.dump(all_queries, outfile)

# def populate_background_json(overall, geog_col):

#     lookup = l[geog_col].upper().replace(" ","_")
#     if lookup not in overall:
#         overall[lookup] = []
    
#     dict_list = overall[lookup]
    
#     lin_found = False
#     for comp_dict in dict_list:
#         if comp_dict['lineage'] == l[lin_col]:
#             comp_dict['count'] += 1
#             lin_found = True
#             break
#     if not lin_found:
#         new_dict = {}
#         new_dict['lineage'] = l[lin_col]
#         new_dict['count'] = 1
#         dict_list.append(new_dict)
        
#     overall[lookup] = dict_list

#     return overall

# def make_background_map_json(config): 

#     lin_col = "lineage" #does this need to be flexible?
#     name_col = config["background_column"]
#     geog_col = config["background_map_column"]

#     start_date, end_date = restrict_by_date(config)

#     overall = defaultdict(list)
#     with open(config["background_csv"]) as f: #this takes about 30 seconds for 1m sequences in Jupyter with no date restrictions. It may require optimisation at some point
#         data = csv.DictReader(f)
#         if config["background_map_date_window"] or config["background_map_date_start"] or config["background_map_date_end"]:
#             for l in data:
#                 if l[geog_col] != "" and l[lin_col] != "":
#                     date = dt.datetime.strptime(l['sample_date'], '%Y-%m-%d').date()
#                     if date >= start_date and date <= end_date:
#                         overall = populate_background_json(overall, geog_col)
#         else:
#             for l in data:
#                 if l[geog_col] != "" and l[lin_col] != "":
#                     overall = populate_background_json(overall, geog_col)
                

#     with open(os.path.join(config["tempdir"],'background_map_data.json'), 'w') as outfile:
#         json.dump(overall, outfile)
