import sys
import os
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
from collections import Counter
import json
import datetime as dt
import requests
import time

def parse_map_file_arg(arg_name, map_file_arg, config): 
    """
    parses map group arguments:
    --background-map-file (default=UK_map if civet_mode == CLIMB, relevant level of world_map (ie adm0 or adm1) if not)
    --query-map-file
    """
    misc.add_arg_to_config(arg_name, map_file_arg, config)

    if config[arg_name]:
        try:
            response = requests.get(config[arg_name])
        except requests.exceptions.MissingSchema:
            print(cyan(f"{config[arg_name]} is an invalid URL.\n"))
            sys.exit(-1)
        if response.status_code != 200:
            sys.stderr.write(cyan(f"Unable to access online resource for {config[arg_name]}. Please ensure you are connected to the internet and the path is valid.\n"))
            sys.exit(-1)
        else:
            map_string = response.text 

        if not config[arg_name].endswith(".json") and not config[arg_name].endswith(".geojson"):
            sys.stderr.write(cyan(f"{config[arg_name]} must be in the format of a geojson. You can use mapshaper.org to convert between file formats.\n"))
            sys.exit(-1)
    else:
        map_string = None

    return map_string

def parse_query_map(query_map_file, longitude_column, latitude_column, found_in_background_data, config):
    """
    parses map group arguments:
    --longitude_column (default=longitude)
    --latitude_column (default=latitude)
    """

    misc.add_arg_to_config("longitude_column",longitude_column,config)
    misc.add_arg_to_config("latitude_column",latitude_column, config)

    if "input_metadata" in config:
        with open(config["input_metadata"]) as f:
            reader = csv.DictReader(f)
            input_fieldnames = reader.fieldnames
    else:
        input_fieldnames = []

    with open(config["background_metadata"]) as f:
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
            with open(config["input_metadata"]) as f:
                data = csv.DictReader(f)
                for line in data:
                    if line[config["input_id_column"]] not in seq_to_longlat_check:
                        seq_to_longlat_check[line[config["input_id_column"]]] = [0,0]
                    if long_input:
                        if line[config['longitude_column']] != "":
                            seq_to_longlat_check[line[config["input_id_column"]]][0] = 1
                    if lat_input:
                        if line[config["latitude_column"]] != "":
                            seq_to_longlat_check[line[config["input_id_column"]]][1] = 1

        data_present = False
        for seq, data in seq_to_longlat_check.items():
            if sum(data) == 2:
                data_present = True
                break
        if not data_present:
            sys.stderr.write(cyan(f"Error: no query with longitude and latitude information contained in the columns provided ({config['longitude_column'], config['latitude_column']}), so map cannot be produced.\n"))
            sys.exit(-1)

        if config["query_map_file"]:
            parse_map_file_arg("query_map_file",query_map_file, config) #don't need to QC this because it doesn't matter what's in it if only the query map
        else:
            if config["civet_mode"] == "CLIMB":
                map_file = "https://viralverity.github.io/civet_geo/uk_map.json"
            else:
                map_file = "https://viralverity.github.io/civet_geo/adm0_global.json"
            
            config["query_map_file"] = map_file


def parse_background_map_options(background_map_file, centroid_file, background_map_date_range, background_map_column, background_map_location, found_in_background_metadata, config):

    """
    parses map group arguments:
    --background-map-date-range (default=no restriction)
    --background-map-location (default=$LOCATION, then adm1 if present, and aggregated_adm2 if civet_mode == CLIMB)
    """

    parse_background_map_column(background_map_column, config)
    parse_date_range(background_map_date_range, config)    
    
    qc_map_file_for_background_map(background_map_file, centroid_file, config)
    check_locations(background_map_location, config)

    qc_centroid_file(config)

    if config["verbose"]:
        print(green("Using the following for mapping background:"))
        print(config["background_map_file"])
        print(config["centroid_file"])

def parse_background_map_column(background_map_column, config):

    misc.add_arg_to_config("background_map_column", background_map_column, config)  

    with open(config["background_metadata"]) as f:
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
        elif "adm0" in background_fieldnames:
            config["background_map_column"] = "adm0"
        elif config["background_location_column"]:
            config["background_map_column"] = config["background_location_column"]
        else:
            sys.stderr.write(cyan(f"Error: no field found in background metadata file for mapping background diversity. Please provide one with -bmcol/--background-map-column.\n") + "\n")
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

            with open(config["background_metadata"]) as f:
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
                    sys.stderr.write(cyan(f"Error: fewer than 50 sequences in column {config['background_map_column']} between dates {start_date} and {end_date}. If you want to map background diversity, consider using a different geographical scale using '-bmloc/--background-map-location' or a larger time period.\n"))
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
    if "input_metadata" in config and config["input_date_column"]:
        with open(config["input_metadata"]) as f:
            data = csv.DictReader(f)
            for l in data:
                if l[config["input_date_column"]] != "": 
                    dates.add(dt.datetime.strptime(l[config["input_date_column"]], '%Y-%m-%d').date())

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

def qc_map_file_for_background_map(background_map_file, centroid_file,config):

    if config["verbose"]:
        print("Beginning checks for background map")

    map_string = parse_map_file_arg("background_map_file", background_map_file, config)
    misc.add_arg_to_config("centroid_file", centroid_file, config)  #needs to happen here so we can check it exists if they have provided a custom map file

    if config["background_map_file"]:

        if not config["centroid_file"]:
            sys.stderr.write(cyan("You have provided a custom background map file, but not a csv containing centroids matching this file. Please provide a csv with the column headers location, longitude and latitude.\n"))
            sys.exit(-1)

        geodata = json.loads(map_string)
        
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
                sys.stderr.write(cyan(f'{config["background_map_column"]} not in default UK map file.\n Options allowed are "suggested_adm2_grouping","adm1" or "adm2".  Alternatively, please provide a custom geojson containing this column to use it using -bmfile/--background-map-file\n'))
                sys.exit(-1)
        else: 
            if config["background_map_column"] == "adm1":
                map_file = "https://viralverity.github.io/civet_geo/adm1_global.json"
            elif config["background_map_column"] == "country" or config["background_map_column"] == "adm0":
                map_file = "https://viralverity.github.io/civet_geo/adm0_global.json"
            else:
                sys.stderr.write(cyan(f"{config['background_map_column']} not in default map file. Please use country/adm0 or adm1 or provide your own shape file using -bmfile/--background-map-file\n"))
                sys.exit(-1)


    config["background_map_file"] = map_file


def check_locations(background_map_location,config):
    #the background files need tidying so that they're all upper case and then we can just check upper case and make them upper case

    acceptable_locations = get_acceptable_locations(config["background_map_file"], config)

    misc.add_arg_to_config("background_map_location", background_map_location,config)

    if config["background_map_location"]:
        lst = config["background_map_location"].split(",")
        new_lst = []
        for i in lst:
            if config["civet_mode"] == "CLIMB":
                if i.upper() not in acceptable_locations:
                    sys.stderr.write(cyan(f'{i} not found in list of acceptable locations to map background lineage diversity.\n Please ensure it is spelt correctly and contains underscores instead of spaces. Please also ensure you are using the correct level of geography by using "--background-map-column". If you still cannot find it and are using default map files, please file a github issue and we will get to it as soon as we can.\n'))
                    sys.exit(-1)
                else:
                    new_lst.append(i.upper())
            else:
                if i.title() not in acceptable_locations:
                    sys.stderr.write(cyan(f'{i} not found in list of acceptable locations to map background lineage diversity.\n Please ensure it is spelt correctly and contains underscores instead of spaces. Please also ensure you are using the correct level of geography by using "--background-map-column". If you still cannot find it and are using default map files, please file a github issue and we will get to it as soon as we can.\n'))
                    sys.exit(-1)
                else:
                    new_lst.append(i.title())

        config["background_map_location"] = new_lst

    else:
        print(green(f"No locations specified for background lineages, so all valid locations in {config['background_map_column']} will be summarised."))
        check_set = set()
        with open(config["background_metadata"]) as f:
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

        missing_locations = set()
        for loc in check_set: 
            if config["civet_mode"] == "CLIMB":
                if loc.upper() not in acceptable_locations:
                    sys.stderr.write(cyan(f'WARNING: {loc} is an invalid location. It will be left out of mapping background diversity\n'))
                    missing_locations.add(loc)
            else:
                if loc.title() not in acceptable_locations:
                    sys.stderr.write(cyan(f'WARNING: {loc} is an invalid location. It will be left out of mapping background diversity\n'))
                    missing_locations.add(loc)

        final_list = []
        for i in check_set:   
            if i not in missing_locations:
                if config["civet_mode"] == "CLIMB":
                    final_list.append(i.upper())
                else:
                    final_list.append(i.title())
        
        config["background_map_location"] = final_list


    if config["verbose"]:
        print("Finished with checks for background map")

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


def qc_centroid_file(config):

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

        for loc in config["background_map_location"]:
            if loc not in centroid_locs:
                sys.stderr.write(cyan(f"{loc}  does not have a centroid associated with it for mapping background diversity. Please add this centroid to the centroid file.\n"))
                sys.exit(-1)

    else:
        if "adm0" in config["background_map_file"]:
            centroid_file = config["adm0_centroids"]
        elif "adm1" in config["background_map_file"]:
            centroid_file = config["adm1_centroids"]
        else:
            centroid_file = config["uk_centroids"]
        
        config["centroid_file"] = centroid_file



### Functions called in report.py ###

def make_query_map_json(config): 
#colour by is dynamic (or will be)
    lat_col = config["latitude_column"]
    long_col = config["longitude_column"]
    name_col = config["background_id_column"] 
 
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

    json_data = json.dumps(all_queries)

    config["start_centre_lat"] = float(start_centre_lat)
    config["start_centre_long"] = float(start_centre_long)

    return json_data


def get_location_information(config, data_for_report):


    with open(config["centroid_file"]) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        
        if "arc_max" in fieldnames:
            get_other_data = True 
        else:
            get_other_data = False
        
        for location in config["background_map_location"]:
            data_for_report[location] = {}

        for row in reader:
            for location in config["background_map_location"]:
                if row['location'] == location:
                    data_for_report[location]["centroids"] = (row["longitude"], row["latitude"])
                    if get_other_data:
                        data_for_report[location]["arc_max"] = float(row["arc_max"])
                        data_for_report[location]["inner_arc_max"] = float(row["inner_arc_max"])
                        data_for_report[location]["text_max"] = float(row["text_max"])
                        data_for_report[location]["start_arc"] = float(row["start_arc"])
                        data_for_report[location]["start_inner_arc"] = float(row["start_inner_arc"])
                        data_for_report[location]["start_text"] = float(row["start_text"])
                    else:
                        data_for_report[location]["arc_max"] = 150
                        data_for_report[location]["inner_arc_max"] = 30
                        data_for_report[location]["text_max"] = 10
                        data_for_report[location]["start_arc"] = 15
                        data_for_report[location]["start_inner_arc"] = 5
                        data_for_report[location]["start_text"] = 5

                    break

    return data_for_report

                

def get_top_ten(counter):

    summary = {}
    top = counter.most_common(10)
    total = sum(list(counter.values()))
    remainder = total
    for lin in top:
        pcent = int(100*(lin[1]/total))
        if pcent >= 5: #otherwise it doesn't get a nice slice of the pie
            remainder-= lin[1]
            summary[lin[0]] = lin[1]
    summary["other"] = remainder
    
    return summary

def make_background_map_json(config): 

    lin_col = "lineage" #does this need to be flexible?
    geog_col = config["background_map_column"]
    wanted_list = set(config["background_map_location"])
    start_date = dt.datetime.strptime(config["start_date"],'%Y-%m-%d').date()
    end_date = dt.datetime.strptime(config["end_date"],'%Y-%m-%d').date()

    locations_all_lins = defaultdict(list)
    with open(config["background_metadata"]) as f: 
        reader = csv.DictReader(f)
        for row in reader:
            if row[geog_col] != "" and row[lin_col] != "":
                if row[geog_col] in wanted_list:
                    date = dt.datetime.strptime(row['sample_date'], '%Y-%m-%d').date()
                    if date >= start_date and date <= end_date:
                        locations_all_lins[row[geog_col]].append(row[lin_col])
                          
    top_ten = defaultdict(dict)
    
    for location, lin_list in locations_all_lins.items():
        counts = Counter(lin_list)
        top_ten[location] = get_top_ten(counts)
       
    overall = []
    
    for location,lineage_counts in top_ten.items():
        for lin, count in lineage_counts.items():
            new_dict = {}
            new_dict["location"] = location
            new_dict["lineage"] = lin
            new_dict["count"] = count
            
            overall.append(new_dict)
    
    json_data = json.dumps(overall)

    return json_data

