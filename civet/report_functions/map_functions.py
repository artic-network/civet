import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
import json
import datetime as dt

def parse_query_map(longitude_column, latitude_column, found_in_background_data, config):
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

def parse_background_map_options(map_file, background_map_date_restriction, background_map_location, found_in_background_metadata, config):

    """
    parses map group arguments:
    --background-map-date-restriction (default=no restriction)
    --background-map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode == CLIMB)
    # --map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode==CLIMB)
    """

    misc.add_arg_to_config("background_map_date_restriction", background_map_date_restriction, config)   
    misc.add_arg_to_config("background_map_location", background_map_location, config)  
    # misc.add_arg_to_config("map_location", map_location, config)  

    with open(config["background_csv"]) as f:
        reader = csv.DictReader(f)
        background_fieldnames = reader.fieldnames

    
    if config["background_map_location"]:
        if config["background_map_location"] not in background_fieldnames: 
            sys.stderr.write(cyan(f"Error: {config['background_map_location']} not found in background metadata file for mapping background diversity\n") + "\n")
            sys.exit(-1)
    
    else:
        if config["civet_mode"] == "CLIMB" and "suggested_adm2_grouping" in background_fieldnames: 
            config["background_map_location"] = "suggested_adm2_grouping"
        elif "adm1" in background_fieldnames:
            config["background_map_location"] = "adm1"
        elif config["location_column"]:
            config["background_map_location"] = config["location_column"]
        else:
            sys.stderr.write(cyan(f"Error: no field found in background metadata file for mapping background diversity. Please provide one with -maploc/--background-map-location.\n") + "\n")
            sys.exit(-1)

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
                        if date >= start_date and date <= end_date and line[config["background_map_location"]] != "":
                            data_count += 1
                        if data_count > 50:
                            enough_data = True
                            break
                if not enough_data:
                    sys.stderr.write(cyan(f"Error: fewer than 50 sequences in column {config['background_map_location']} between dates {start_date} and {end_date}. If you want to map background diversity, consider using a different geographical scale using '-maploc/--background-map-location' or a larger time period.\n"))
                    sys.exit(-1)
        else:
            sys.stderr.write(cyan(f"Error: Date restriction defined for background lineage diversity mapping, but no date column provided. Please use '-bdate/--background-date-column' to specify this column. \n") + "\n")
            sys.exit(-1)
    else:
        print(green("No date restriction provided, analysing background lineage diversity for the whole pandemic.\n"))
        config["start_date"] = dt.datetime(2019,12,1).date()
        config["end_date"] = dt.datetime.today().date()

    # check_shapefile(map_file,config)


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

def parse_map_file(map_file, config):
    #winding in the geojsons

    #still need to make UK shapefile - no, this exists in b117 paper
    #work out what's going on with the LFS stuff with the larger shapefiles - maybe just add them to a website?

    """
    parses map group arguments:
    --map-file (default=UK_map if civet_mode == CLIMB, relevant level of world_map (ie adm0 or adm1) if not)
    """

    misc.add_arg_to_config("map_file", map_file, config)

    if config["map_file"]:
        if not os.path.exists(config["map_file"]):
            sys.stderr.write(cyan(f"{config['map_file']} cannot be found. Please check the path and try again.\n"))
            sys.exit(-1)
        if not config["map_file"].endswith(".json") and not config['map_file'].endswith(".geojson"):
            sys.stderr.write(cyan(f"{config['map_file']} must be in the format of a geojson. You can use mapshaper.org to convert between file formats.\n"))
            sys.exit(-1)

        with open(config["map_file"]) as f:
            geodata = json.load(f)
        
        headers = geodata["features"][0]["properties"].keys()
        
        if config["background_map_location"] not in headers:
            sys.stderr.write(cyan(f"{config['background_map_location']} not found in custom shapefile.\n"))
            sys.exit(-1)
        else:
            acceptable_locations = []
            for item in geodata['features']:
                acceptable_locations.append(item["properties"][config["background_map_location"]])
    
    else:
        if config["civet_mode"] == "CLIMB":
            uk_cols = ["suggested_adm2_grouping", "adm1", "adm2"]
            map_file = config["uk_map_file"]
            if config["background_map_location"] not in uk_cols:
                sys.stderr.write(cyan(f'{config["background_map_location"]} not in default UK shapefile. Please provide a custom geojson containing this column to use it using --map-file\n'))
                sys.exit(-1)
        else:
            if config["background_map_location"] == "adm1":
                map_file = config["adm1_global_file"]
            elif config["background_map_location"] == "country" or config["background_map_location"] == "adm0" or config["background_map_location"] == "ISO":
                map_file = config["adm0_global_file"]
            else:
                sys.stderr.write(cyan(f"{config['background_map_location']} not in default shapefiles. Please use country/adm0 or adm1 or provide your own shape file using --map-file\n"))
                sys.exit(-1)

        acceptable_locations = get_acceptable_locations(map_file, config)

    check_set = set()
    with open(config["background_csv"]) as f:
        data = csv.DictReader(f)
        for line in data:
            if config["background_date_column"]:
                date_value = line[config["background_date_column"]]
                if date_value != "":
                    date = dt.datetime.strptime(line[config["background_date_column"]]).date()
                    if date >= config["start_date"] and date <= config["end_date"]:
                        check_set.add(line[config["background_map_location"]])
            else:
                check_set.add(line[config["background_map_location"]])

    for loc in check_set:
        if loc not in acceptable_locations:
            if config["map_file"]:
                sys.stderr(f'{loc} is an invalid location. Please ensure that the metadata values match up to the shapefile you have provided.\n')
                sys.exit(-1)
            elif config['civet_mode'] == "CLIMB":
                sys.stderr(f'{loc} is an invalid location. If you are using the default background metadata, please contact Verity Hill (verity.hill@ed.ac.uk)')
                sys.exit(-1)
            else:
                sys.stderr(f"{loc} isn't in our shapefile. Please see a list of currently accepted locations here: [link]. If you can't find your country's data on that list, please open a github issue on the civet repo aand we will get to it as soon as we can.\n")
                sys.exit(-1)

    config["map_file"] = map_file

    return

def get_acceptable_locations(map_file, config):

    if config["civet_mode"] == "CLIMB":
        if config["background_map_location"] != "adm1":
            acceptable_locations = ["Scotland", "Wales", "Northern Ireland", "England"]
        else:
            acceptable_locations = []
            with open(config["uk_acceptable_values"]) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    acceptable_locations.append(reader[config["background_map_location"]])

    else:
        acceptable_locations = set()
        with open(config["global_accepted_values"]) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                acceptable_locations.add(reader[config["background_map_location"]])
    
    return acceptable_locations
                

# def make_query_map_json(config, metadata): 
# #plan would be to just plot this on the world map JSON that you can then zoom in on wherever required
# #Then have interactivity for colour-by. Easy enough to add in here to the JSON and QC checks if that's too complicated to have interactive
# #Possible extension: can generate long/lat 

#     lat_col = config["latitude_column"]
#     long_col = config["longitude_column"]
#     name_col = config["background_column"] 

#     overall_dict = defaultdict(dict)

#     all_queries = {}
#     all_queries['queries'] = []

#     with open(metadata) as f:
#         data = csv.DictReader(f)
#         for l in data:
#             if l["query_boolean"] == "TRUE":
#                 per_seq_dict = {}
#                 per_seq_dict["sequence_name"] = l[name_col]
#                 per_seq_dict["latitude"] = l[lat_col]
#                 per_seq_dict["longitude"] = l[long_col]
                
#                 all_queries['queries'].append(per_seq_dict)

#     with open(os.path.join(config["tempdir"],'query_map_data.json'), 'w') as outfile:
#         json.dump(all_queries, outfile)

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
