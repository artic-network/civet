import sys
import os
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
import json
import datetime as dt


def parse_map_file_arg(map_file_arg, arg_name, config): 
#this needs to change because it has to be an online resource
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
                map_file = "https://viralverity.github.io/civet_geo/uk_map.json"
            else:
                map_file = "https://viralverity.github.io/civet_geo/adm0_global.json"
            
            config["query_map_file"] = map_file


def parse_background_map_options(background_map_file, background_map_date_range, background_map_column, background_map_location, found_in_background_metadata, config):

    """
    parses map group arguments:
    --background-map-date-range (default=no restriction)
    --background-map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode == CLIMB)
    """

    parse_background_map_column(background_map_column, config)
    parse_date_range(background_map_date_range, config)    
    
    map_file, acceptable_locations = qc_map_file_for_background_map(background_map_file, config)

    config["background_map_file"] = map_file

    qc_centroid_file(config, acceptable_locations)
    qc_background_map_location(background_map_location,acceptable_locations, config)

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

def parse_date_range(background_map_date_range, config):

    misc.add_arg_to_config("background_map_date_range", background_map_date_range, config)   

    if config["background_map_date_range"]:
        if config["background_date_column"]:
            date_range = config["background_map_date_range"].split(":")
            if len(date_range) > 1:
                try:
                    start_date = dt.datetime.strptime(date_range[0], "%Y-%m-%d").date()
                except:
                    sys.stderr.write(cyan("Start date in background map date range in incorrect format. Please use YYYY-MM-DD"))
                    sys.exit(-1)
                try:
                    end_date = dt.datetime.strptime(date_range[1], "%Y-%m-%d").date()
                except:
                    sys.stderr.write(cyan("End date in background map date range in incorrect format. Please use YYYY-MM-DD"))
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

def qc_map_file_for_background_map(config):

    if config["verbose"]:
        print("Beginning checks for background map")

    parse_map_file_arg(background_map_file, "background_map_file", config)
    misc.add_arg_to_config("centroid_file", centroid_file, config)  #needs to happen here so we can check it exists if they have provided a custom map file

    if config["background_map_file"]:

        if not config["centroid_file"]:
            sys.stderr.write(cyan("You have provided a custom background map file, but not a csv containing centroids matching this file. Please provide a csv with the column headers location, longitude and latitude.\n"))
            sys.exit(-1)

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
            map_file = "https://viralverity.github.io/civet_geo/uk_map.json"
            uk_cols = ["suggested_adm2_grouping", "adm1", "adm2"]
            if config["background_map_column"] not in uk_cols:
                sys.stderr.write(cyan(f'{config["background_map_column"]} not in default UK map file.\n Options allowed are "suggested_adm2_grouping","adm1" or "adm2".  Alternatively, please provide a custom geojson containing this column to use it using --background-map-file\n'))
                sys.exit(-1)
        else: 
            if config["background_map_column"] == "adm1":
                map_file = "https://viralverity.github.io/civet_geo/adm1_global.json"
            elif config["background_map_column"] == "country" or config["background_map_column"] == "adm0" or config["background_map_column"] == "ISO":
                map_file = "https://viralverity.github.io/civet_geo/adm0_global.json"
            else:
                sys.stderr.write(cyan(f"{config['background_map_column']} not in default map file. Please use country/adm0 or adm1 or provide your own shape file using --background-map-file\n"))
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

    for loc in check_set: 
        if loc not in acceptable_locations:
            if config["background_map_file"]:
                sys.stderr.write(cyan(f'{loc} is an invalid location. Please ensure that the metadata values match up to the map file you have provided.\n'))
                sys.exit(-1)
            elif config['civet_mode'] == "CLIMB":
                sys.stderr.write(cyan(f'{loc} is an invalid location. If you are using the default background metadata, please contact Verity Hill (verity.hill@ed.ac.uk)\n'))
                sys.exit(-1)
            else:
                sys.stderr.write(cyan(f"{loc} isn't in our map file. Please see a list of currently accepted locations here: [link]. If you can't find your country's data on that list, please open a github issue on the civet repo and we will get to it as soon as we can.\n"))
                sys.exit(-1)

    if config["verbose"]:
        print("Finished with checks for background map")

    return map_file, acceptable_locations

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


def qc_centroid_file(config, acceptable_locations):

    if config["centroid_file"]:
        if not os.path.exists(config["centroid_file"]):
            sys.stderr.write(f"Cannot find centroid file at {config['centroid_file']}.\n")
            sys.exit(-1)
        required_headers = ["location", "latitude", "longitude"]
        with open(config["centroid_file"]) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            for i in required_headers:
                if i not in fieldnames:
                    sys.stderr.write(f"{i} is not in the provided centroid file. Please provide a centroid file as a csv with the headers location, latitude and longitude.\n")
                    sys.exit(-1)
            centroid_locs = []
            for row in reader:
                centroid_locs.append(row["location"])

        for loc in acceptable_locations:
            if loc not in centroid_locs:
                sys.stderr.write(cyan(f"{loc} is in the provided geojson but does not have a centroid associated with it. Please add this centroid to the centroid file.\n"))
                sys.exit(-1)

    else:
        if "adm0" in config["background_map_file"]:
            centroid_file = config["adm0_centroids"]
        elif "adm1" in config["background_map_file"]:
            centroid_file = config["adm1_centroids"]
        else:
            centroid_file = config["uk_centroids"]
        
        config["centroid_file"] = centroid_file


def qc_background_map_location(background_map_location,acceptable_locations, config):

    misc.add_arg_to_config(background_map_location, "background_map_location", config)

    if config["background_map_location"]:
        lst = config["background_map_location"].split(",")
        for i in lst:
            if i not in acceptable_locations:
                sys.stderr.write(cyan(f'{i} not found in list of acceptable locations to map background lineage diversity.\n'))
                sys.exit(-1)
        config["background_map_location"] = lst
    else:
        config["background_map_location"] = acceptable_locations

### Functions called in report.py ###

def make_query_map_json(config): 
#colour by is dynamic (or will be)
    lat_col = config["latitude_column"]
    long_col = config["longitude_column"]
    name_col = config["background_column"] 
 
    all_queries = []

    with open(config["query_metadata"]) as f:
        data = csv.DictReader(f)
        for l in data:
            per_seq_dict = {}
            if l[lat_col] != "" and l[long_col] != "":
                per_seq_dict["sequence_name"] = l[name_col]
                per_seq_dict["latitude"] = l[lat_col]
                per_seq_dict["longitude"] = l[long_col]
                start_centre_lat = l[lat_col]
                start_centre_long = l[long_col]
                
                all_queries.append(per_seq_dict)

    json_name = os.path.join(config["tempdir"], 'query_map_data.json')

    with open(json_name, 'w') as outfile:
        json.dump(all_queries, outfile)

    config["start_centre_lat"] = float(start_centre_lat)
    config["start_centre_long"] = float(start_centre_long)

    return json_name


def get_centroids(config):

    centroid_dict = {}
    with open(config["centroid_file"]) as f:
        reader = csv.DictReader(f)
        for row in reader:
            centroid_dict[row["location"]] = (row["longitude"], row["latitude"])

    return centroid_dict

                

def get_top_ten(counter):

    summary = {}
    top = counter.most_common(10)
    total = sum(list(counter.values()))
    print(top)
    remainder = total
    for lin in top:
        pcent = int(100*(lin[1]/total))
        if pcent >= 1:
            remainder-= lin[1]
            summary[lin[0]] = lin[1]
    summary["other"] = remainder
    
    return summary

def make_background_map_json(config): 

    lin_col = "lineage" #does this need to be flexible?
    geog_col = config["background_map_column"]
    wanted_list = set(config["background_map_location"])

    locations_all_lins = defaultdict(list)
    with open(config["background_csv"]) as f: 
        data = csv.DictReader(f)
        for l in tqdm.tqdm(data):
            if l[geog_col] != "" and l[lin_col] != "":
                if l[geog_col] in wanted_list:
                    date = dt.datetime.strptime(l['sample_date'], '%Y-%m-%d').date()
                    if date >= config["start_date"] and date <= config["end_date"]:
                        locations_all_lins[l[geog_col]].append(l[lin_col])
                        
                        
    top_ten = defaultdict(dict)
    
    for location, lin_list in locations_all_lins.items():
        counts = Counter(lin_list)
        top_ten[location] = get_top_ten(counts)
        
       
    overall = []
    
    for location,lineage_counts in location_top_tens.items():
        for lin, count in lineage_counts.items():
            new_dict = {}
            new_dict["location"] = location
            new_dict["lineage"] = lin
            new_dict["count"] = count
            
            overall.append(new_dict)
    
    with open('background_map_data.json', 'w') as outfile:
        json.dump(overall, outfile)

