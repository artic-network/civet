import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt
from collections import defaultdict
import json

## Should just take any map file - should have a check that some of the sequences in match to 


##Think about what the defaults are - for local lins, could be adm1. 
# If in the UK, would want it to be the "suggested_aggregated" or w/e that col is called. Could have it as adm1 and tell people using it we advise agg_adm2



def parse_map_options(metadata, map_queries, map_background, latitude_column, longitude_column, background_map_column,background_map_date_window, background_map_date_start,background_map_date_end, background_map_date_column, config):

    misc.add_arg_to_config("map_background", map_background, config)
    misc.add_arg_to_config("background_map_column", background_map_column, config)
    misc.add_arg_to_config("background_map_date_start", background_map_date_start, config)
    misc.add_arg_to_config("background_map_date_end", background_map_date_end, config)
    misc.add_arg_to_config("background_map_date_window", background_map_date_window, config)
    misc.add_arg_to_config("background_map_date_column", background_map_date_column, config)

    misc.add_arg_to_config("map_queries",map_queries,config)
    misc.add_arg_to_config("longitude_column",longitude_column,config)
    misc.add_arg_to_config("latitude_column", latitude_column, config)

    if config["map_background"] or config["map_queries"]:
        with open(metadata) as f:
            data = csv.DictReader(f)
            headers = data.fieldnames
            if config["map_background"]:
                if not config["background_map_column"]:
                    sys.stderr.write(cyan(f'Error: --map-background provided, but no column containing the geographical data to summarise background lineages by provided.\n') + "Please specify using '--background-map-column'.\n")
                    sys.exit(-1)
                elif config["background_map_column"] not in headers:
                    sys.stderr.write(cyan(f"Error: {config['background_map_column']} not found in metadata file for mapping background diversity\n") + "\n")
                    sys.exit(-1)
                else:

                    if config["background_map_date_start"]:
                        try:
                            dt.datetime.strptime(config["background_map_date_start"], "%Y-%m-%d").date()
                        except:
                            sys.stderr.write(cyan(f"Error: Please provide 'background_map_date_start' argument in YYYY-MM-DD format.\n"))
                    if config["background_map_date_end"]:
                        try:
                            dt.datetime.strptime(config["background_map_date_end"], "%Y-%m-%d").date()
                        except:
                            sys.stderr.write(cyan(f"Error: Please provide 'background_map_date_end' argument in YYYY-MM-DD format.\n"))


                    if config["background_map_date_start"] or config["background_map_date_end"] or config["background_map_date_window"]:
                        if not config["background_map_date_column"]:
                            sys.stderr.write(cyan(f"Error: no date column provided to restrict background lineage diversity analysis.\n")+ "Please either provide a date column or remove the date restriction arguments to show background lineage diversity over the whole pandemic. \n")
                            sys.exit(-1)
                        elif config["background_map_date_column"] not in data.fieldnames:
                            sys.stderr.write(cyan(f"Error: {config['background_map_date_column']} not found in metadata file for restricting the mapping of background diversity to a specific time window\n"))
                            sys.exit(-1)
                        else:
                            line_count = 0
                            date_present = False
                            for line in data:
                                line_count += 1
                                if line[config["background_map_date_column"]] != "":
                                    misc.check_date_format(line[config["background_map_date_column"]], line_count, config["background_map_date_column"])
                                    date_present = True
                            if not date_present:
                                sys.stderr.write(cyan(f"Error: no date information in {config["background_map_date_column"]}\n") + "Please provide a column with date information, or remove date restriction arguments to show background lineage diversity over the whole pandemic. \n")
                    
                    else:
                        print(green("No date restriction provided, analysing background lineage diversity for the whole pandemic.\n"))

                    check_shapefile(config, "background", metadata, shapefile_location)

            if config["map_queries"]:
                if not config["longitude_column"] or not config["latitude_column"]:
                    sys.stderr.write(cyan(f'Error: --map-queries provided, but no columns containing longitude and latitude coordinates specified.\n') + "Please specify using '--latitude-column' and '--longitude-column'\n")
                    sys.exit(-1)
                elif config["longitude_column"] not in headers:
                    sys.stderr.write(cyan(f"Error: {config['longitude_column']} not found in metadata file for mapping queries\n") + "\n")
                    sys.exit(-1)
                elif config["latitude_column"] not in headers:
                    sys.stderr.write(cyan(f"Error: {config['latitude_column']} not found in metadata file for mapping queries\n") + "\n")
                    sys.exit(-1)
                
                check_shapefile(config, "queries", metadata, shapefile_location)


    if (config["background_map_date_column"] or config["background_map_column"] or config["background_map_date_start"] or config["background_map_date_end"] or config["background_map_date_window"]) and not config["map_background"]:
        sys.stderr.write(cyan(f'Arguments for mapping background lineage diversity, but --map-background not used.\n') + 'Please also use --map-background to use this function\n')
    if config["query_map_column"] and not config["map_queries"]:
        sys.stderr.write(cyan(f'Arguments for mapping queries diversity, but --map-queries not used.\n') + 'Please also use --map-queries to use this function\n')

   return config



def make_query_map_json(config, metadata):

    lat_col = config["latitude_column"]
    long_col = config["longitude_column"]
    name_col = config["data_column"] 

    overall_dict = defaultdict(dict)

    all_queries = {}
    all_queries['queries'] = []

    with open(metadata) as f:
        data = csv.DictReader(f)
        for l in data:
            if l["query_boolean"] == "TRUE":
                per_seq_dict = {}
                per_seq_dict["sequence_name"] = l[name_col]
                per_seq_dict["latitude"] = l[lat_col]
                per_seq_dict["longitude"] = l[long_col]
                
                all_queries['queries'].append(per_seq_dict)

    with open(os.path.join(config["tempdir"],'query_map_data.json', 'w') as outfile:
        json.dump(all_queries, outfile)


def make_background_map_json(config, background_metadata): #this needs the full background metadata, not the query metadata




def restrict_by_date(config, metadata):





def check_shapefile(config, map_type, metadata, utils_dir):
## How are we doing shape files - this is probably the biggest unsolved problem for now
##Maybe:
    #if col is country, then can use one that we have
    #could we store adm1 for every country? Look into that. Will have to do a bit of optimising there probably 
    #if col is adm2 and UK, no problem, we'll provide that. 
    #if col is adm2 and not UK, they have to provide a shapefile
    #will also have to think about the title in the shapefile

#Docs for this are going to need to be good!




