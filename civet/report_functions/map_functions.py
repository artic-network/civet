import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
import json


##Think about what the defaults are - for local lins, could be adm1. 
# If in the UK, would want it to be the "suggested_aggregated" or w/e that col is called. Could have it as adm1 and tell people using it we advise agg_adm2

#will have to have levels of defaults for location - agg_adm2 if present, adm1 if present, country if present. If none present and not given, then crap out


def parse_map_file():

    return

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

# def parse_map_options(map_file, longitude_column, latitude_column, background_map_date_restriction, config):

#     misc.add_arg_to_config("background_map_date_restriction", background_map_date_restriction, config)    


#     if config["map_background"]:

#         if not config["background_map_column"]:
#             sys.stderr.write(cyan(f'Error: --map-background provided, but no column containing the geographical data to summarise background lineages by provided.\n') + "Please specify using '--background-map-column'.\n")
#             sys.exit(-1)
#         # elif config["background_map_column"] not in headers: #possibly don't need this check - depends on how we zoom in. If it's all interactive, we don't need it
#         #     sys.stderr.write(cyan(f"Error: {config['background_map_column']} not found in metadata file for mapping background diversity\n") + "\n")
#         #     sys.exit(-1)
#         else:

#             if config["background_map_date_start"]:
#                 try:
#                     dt.datetime.strptime(config["background_map_date_start"], "%Y-%m-%d").date()
#                 except:
#                     sys.stderr.write(cyan(f"Error: Please provide 'background_map_date_start' argument in YYYY-MM-DD format.\n"))
#             if config["background_map_date_end"]:
#                 try:
#                     dt.datetime.strptime(config["background_map_date_end"], "%Y-%m-%d").date()
#                 except:
#                     sys.stderr.write(cyan(f"Error: Please provide 'background_map_date_end' argument in YYYY-MM-DD format.\n"))

#             if config["background_map_date_start"] or config["background_map_date_end"] or config["background_map_date_window"]:
#                 if not config["background_map_date_column"]:
#                     sys.stderr.write(cyan(f"Error: no date column provided to restrict background lineage diversity analysis.\n")+ "Please either provide a date column or remove the date restriction arguments to show background lineage diversity over the whole pandemic. \n")
#                     sys.exit(-1)
#                 else:
#                     with open(metadata) as f:
#                         data = csv.DictReader(f)
#                         headers = data.fieldnames
#                         if config["background_map_date_column"] not in data.fieldnames:
#                             sys.stderr.write(cyan(f"Error: {config['background_map_date_column']} not found in query metadata file for restricting the mapping of background diversity to a specific time window\n"))
#                             sys.exit(-1)
#                         else:
#                             line_count = 0
#                             date_present = False
#                             for line in data:
#                                 line_count += 1
#                                 if line[config["background_map_date_column"]] != "":
#                                     misc.check_date_format(line[config["background_map_date_column"]], line_count, config["background_map_date_column"])
#                                     date_present = True
#                             if not date_present:
#                                 sys.stderr.write(cyan(f"Error: no date information in {config['background_map_date_column']}\n") + "Please provide a column with date information, or remove date restriction arguments to show background lineage diversity over the whole pandemic. \n")
                
#                     with open(config["background_csv"]) as f:
#                         data = csv.DictReader(f)
#                         bg_headers = data.fieldnames()
#                         if config["background_map_column"] not in bg_headers: 
#                             sys.stderr.write(cyan(f"Error: {config['background_map_column']} not found in background metadata file for mapping background diversity\n") + "\n")
#                             sys.exit(-1)
#                         elif "lineage" not in headers:
#                             sys.stderr.write(cyan(f"Error: 'lineage' not found in metadata file for mapping background diversity\n") + "\n")
#                             sys.exit(-1)
#                         elif config["background_map_date_column"] not in bg_headers:
#                             sys.stderr.write(cyan(f"Error: {config['background_map_date_column']} not found in background metadata file for restricting the mapping of background diversity to a specific time window\n"))
#                             sys.exit(-1)

#                         else:
#                             line_count = 0
#                             for l in data:
#                                 line_count += 1
#                                 if line[config["background_map_date_column"]] != "":
#                                     misc.check_date_format(line[config["background_map_date_column"]], line_count, config["background_map_date_column"])
#                                     date_present = True
#                             if not date_present:
#                                 sys.stderr.write(cyan(f"Error: no date information in {config['background_map_date_column']} in background metadata.\n") + "Please provide a column with date information, or remove date restriction arguments to show background lineage diversity over the whole pandemic. \n")


#             else:
#                 print(green("No date restriction provided, analysing background lineage diversity for the whole pandemic.\n"))

#                 check_shapefile(config, "background", metadata, shapefile_location)
        
        
            
#     if config["map_queries"]:
#         with open(metadata) as f:
#             data = csv.DictReader(f)
#             headers = data.fieldnames
#             if not config["longitude_column"] or not config["latitude_column"]:
#                 sys.stderr.write(cyan(f'Error: --map-queries provided, but no columns containing longitude and latitude coordinates specified.\n') + "Please specify using '--latitude-column' and '--longitude-column'\n")
#                 sys.exit(-1)
#             elif config["longitude_column"] not in headers:
#                 sys.stderr.write(cyan(f"Error: {config['longitude_column']} not found in metadata file for mapping queries\n") + "\n")
#                 sys.exit(-1)
#             elif config["latitude_column"] not in headers:
#                 sys.stderr.write(cyan(f"Error: {config['latitude_column']} not found in metadata file for mapping queries\n") + "\n")
#                 sys.exit(-1)
#             else:
#                 data_present = False
#                 for line in data:
#                     if line["query_boolean"] == "TRUE" and line[config["latitude_column"]] != "" and line[config["longitude_column"]] != "":
#                         data_present = True
#                         break
#                 if not data_present:
#                     print(cyan("Error: no query with longitude and latitude information, so map cannot be produced.\n"))

                
#             check_shapefile(config, "queries", metadata, shapefile_location)


#     # if (config["background_map_date_column"] or config["background_map_column"] or config["background_map_date_start"] or config["background_map_date_end"] or config["background_map_date_window"]) and not config["map_background"]:
#     #     sys.stderr.write(cyan(f'Arguments for mapping background lineage diversity, but --map-background not used.\n') + 'Please also use --map-background to use this function\n')
#     # if config["query_map_column"] and not config["map_queries"]:
#     #     sys.stderr.write(cyan(f'Arguments for mapping queries, but --map-queries not used.\n') + 'Please also use --map-queries to use this function\n')
#     #commented out because we want to be able to have defaults for these




# def make_query_map_json(config, metadata): 
# #plan would be to just plot this on the world map JSON that you can then zoom in on wherever required
# #Then have interactivity for colour-by. Easy enough to add in here to the JSON and QC checks if that's too complicated to have interactive
# #Possible extension: can generate long/lat 

#     lat_col = config["latitude_column"]
#     long_col = config["longitude_column"]
#     name_col = config["data_column"] 

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
#     name_col = config["data_column"]
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


# def restrict_by_date(config):

#     start_date = None
#     end_date = None

#     if config["background_map_date_window"]:
#         date_diff = dt.timedelta(config["background_map_date_window"])
#         dates = set()
#         with open(config["query_metadata"]) as f:
#             data = csv.DictReader(f)
#             for l in data:
#                 if l['query_boolean'] == 'TRUE':
#                     if l[config['background_map_date_column'] != "": #default will be sample_date, it's just in case they provide their own metadata
#                         dates.add(dt.datetime.strptime(l['sample_date'], '%Y-%m-%d').date())
                
#         earliest = sorted(dates)[0]
#         latest = sorted(dates)[-1]

#         start_date = earliest - date_diff
#         end_date = latest + date_diff

#         return start_date, end_date

#     if config["background_map_date_start"]:
#         start_date = config["background_map_date_start"]

#     if config["background_map_date_end"]:
#         end_date = config["background_map_date_end"]

#     if config["background_map_date_start"] and not config["background_map_date_end"]:
#         end_date = dt.datetime.today().date()
#     elif config["background_map_date_end"] and not config["background_map_date_start"]:
#         start_date = dt.date(2019,1,12)

#     return start_date, end_date


# def check_shapefile(config, map_type, metadata, utils_dir):
## How are we doing shape files - this is probably the biggest unsolved problem for now
##Maybe:
    #if col is country, then can use one that we have
    #could we store adm1 for every country? Look into that. Will have to do a bit of optimising there probably 
    #if col is adm2 and UK, no problem, we'll provide that. 
    #if col is adm2 and not UK, they have to provide a shapefile
    #will also have to think about the title in the shapefile

#Docs for this are going to need to be good!




