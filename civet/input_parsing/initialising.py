from datetime import date
import os
import sys
import yaml
from civet.utils.log_colours import green,cyan

import civet.utils.custom_logger as custom_logger
import civet.utils.log_handler_handle as lh
from civet.utils import misc

def get_defaults():
    today = date.today()
    default_dict = {
                    "date": today,# date investigation was opened
                    "authors": "", # List of authors, affiliations and contact details
                    "max_queries":5000,

                    # Initialising data variables
                    "num_seqs":0,
                    "datadir": os.getenv('DATADIR'),
                    "background_csv":False,
                    "background_fasta":False,
                    "background_SNPs":False,
                    "background_tree":False,

                    # Output defaults
                    "output_prefix":"civet",
                    "output_data":False,
                    "datestamp":False,
                    "no_temp":False,
                    "overwrite":False,

                    # Search defined by metadata column
                    "from_metadata":False,

                    # Columns to match
                    "input_column":"name",
                    "background_column":"sequence_name",
                    "fasta_column":False,
                    
                    # reference coordinates to pad until and from
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674,

                    # catchment config
                    "catchment_size":100,
                    "distance":False,
                    "distance_up":2,
                    "distance_down":2,
                    "distance_side":2,
                    "collapse_threshold":1,

                    "downsample":["mode=random"],
                    "mode":False,
                    "factor":False,
                    "downsample_field":False,
                    "downsample_column":False,
                    "max_memory": 8,

                    # QC standards for input fasta file
                    "max_ambiguity":0.5,
                    "min_length":20000,

                    #Tree options
                    "tree_annotations": "lineage,country",

                    #Report options
                    "colour_theme": "#7178bc",
                    "colour_map": "#86b0a6,#565b7b,#e9c46a,#e1998a,#d19c2c,#264653,#4888a3,#c77939,#b56576,#eb7e83,#F46D43,#9e2a2b,#84a98c",
                    "report_content":["1,2,3,4,5"],
                    "report_preset":False,
                    "report_title":f"civet report",
                    "anonymise":False,
                    "report_column": False,
                    "timeline_dates":False, #default is date_column or background_date_column if they are present
                    "timeline_colours":False,
                    "table_content":False,
                    "background_date_column":False, #default is sample_date if present in the background csv, then date_column if present
                    "date_column":False, #default is sample_date if present in the input csv
                    "location_column":False, #default is country if present in the background csv

                    #Map options
                    "latitude_column":"latitude",
                    "longitude_column":"longitude",
                    "query_map_file":False, #if UK, will find the geojson containing aggregated_adm2s, if global will find adm0 or adm1 geojson
                    "background_map_file": False, #Same as above
                    "centroid_file":False, #will select appropriate centroid file for provided jsons as default
                    "background_map_column":False, #default is location_column
                    "background_map_location":False, #default is all valid locations in background metadata
                    # "map_column":False,
                    "background_map_date_range":False, #default is no restriction, will do background from 2019-11-01 to today

                    # misc defaults
                    "civet_mode": os.getenv('CIVETMODE'),
                    "threads":1,
                    "verbose":False
                    }
    return default_dict

def check_configfile(cwd,config_arg):
    configfile = os.path.join(cwd,config_arg)

    ending = configfile.split(".")[-1]

    if ending not in ["yaml","yml"]:
        sys.stderr.write(cyan(f'Error: config file {configfile} must be in yaml format.\n'))
        sys.exit(-1)
    
    elif not os.path.isfile(configfile):
        sys.stderr.write(cyan(f'Error: cannot find config file at {configfile}\n'))
        sys.exit(-1)
    else:
        print(green(f"Input config file:") + f" {configfile}")
        return configfile

def arg_dict(config):
    arguments = {
                # igroup args
                "i":"input_csv",
                "input_csv":"input_csv",
                "icol":"input_column",
                "input_column":"input_column",
                "ids":"id_string",
                "id_string":"id_string",
                "fm":"from_metadata",
                "from_metadata":"from_metadata",

                "f":"fasta",
                "fasta":"fasta",
                "n":"max_ambiguity",
                "max_ambiguity":"max_ambiguity",
                "l":"min_length",
                "min_length":"min_length",

                # dgroup args
                "datadir":"datadir",
                "DATADIR":"datadir",
                "d":"datadir",
                "background_csv":"background_csv",
                "bc":"background_csv",
                "background_fasta":"background_fasta",
                "bf":"background_fasta",
                "background_SNPs":"background_SNPs",
                "bSNP":"background_SNPs",
                "background_tree":"background_tree",
                "bt":"background_tree",
                "bcol":"background_column",
                "background_column":"background_column",
                "fcol":"fasta_column",
                "fasta_column":"fasta_column",

                # ogroup args
                "o":"outdir",
                "outdir":"outdir",
                "p":"output_prefix",
                "output_prefix":"output_prefix",
                "datestamp":"datestamp",
                "overwrite":"overwrite",
                "output_data":"output-data",
                "no-temp":"no-temp",
                "tempdir":"tempdir",
                "temp":"tempdir",

                # agroup args
                "ts":"trim_start",
                "te":"trim_end",
                "trim_start":"trim_start",
                "trim_end":"trim_end",
                "cs":"catchment_size",
                "catchment_size":"catchment_size",
                "max_queries":"max_queries",
                "mq":"max_queries",
                "max_memory":"max_memory",
                "mem":"max_memory",
                "downsample":"downsample",
                "ds":"downsample",
                "distance":"distance",
                "dSNP":"distance",
                "distance_up":"distance_up",
                "distance_down":"distance_down",
                "distance_side":"distance_side",

                #rgroup args
                "colour_theme":"colour_theme",
                "ct":"colour_theme",
                "colour_map":"colour_map",
                "cmap":"colour_map",
                "report_title":"report_title",
                "rt":"report_title",
                "report_content":"report_content",
                "rc":"report_content",
                "report_preset":"report_preset",
                "rp":"report_preset",

                "alt":"report_column",
                "report_column":"report_column",

                "anonymise":"anonymise",
                "anonymize":"anonymise",

                "td":"timeline_dates",
                "timeline_dates":"timeline_dates",
                "timeline_colours":"timeline_colours",
                "timeline_colors":"timeline_colours",
                "table_content":"table_content",

                "background_date_column":"background_date_column",
                "bdate":"background_date_column",
                "date_column":"date_column",
                "date":"date_column",

                "location_column":"location_column",
                "lcol":"location_column",

                "tree_annotations":"tree_annotations",

                #mgroup args
                "lat":"latitude_column",
                "latitude_column":"latitude_column",
                "long":"longitude_column",
                "longitude_column":"longitude_column",
                "query_map_file":"query_map_file",
                "qmfile":"query_map_file",
                
                "background_map_date_range":"background_map_date_range",
                "bmrange":"background_map_date_range",
                "background_map_column":"background_map_column",
                "bmcol":"background_map_column",
                "background_map_file":"background_map_file",
                "background_map_location":"background_map_location",
                "bmloc":"background_map_location",
                "bmfile":"background_map_file",
                "centroid_file":"centroid_file",

                # misc group args
                "civet_mode":"civet_mode",
                "t":"threads",
                "threads":"threads",
                "v":"verbose",
                "verbose":"verbose"
    }
    for i in config:
        arguments[i] = i
    return arguments

def load_yaml(f):
    try:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    except:
        sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
        sys.exit(-1)
    return input_config

def return_path_keys():
    return ["input_csv","fasta","background_csv","background_fasta","background_tree","background_SNPs","datadir","outdir","tempdir"]

def setup_absolute_paths(path_to_file,value):
    return os.path.join(path_to_file,value)


def parse_yaml_file(configfile,configdict):
    overwriting = 0
    path_keys = return_path_keys()

    path_to_file = os.path.abspath(os.path.dirname(configfile))
    configdict["input_path"] = path_to_file
    valid_keys = arg_dict(configdict)

    invalid_keys = []
    with open(configfile,"r") as f:
        input_config = load_yaml(f) # try load file else exit with msg

        for key in input_config:
            value = input_config[key]
            if value == None: # dont count blank entries
                pass
            else:
                clean_key = key.lstrip("-").replace("-","_").rstrip(" ").lstrip(" ").lower()

                if clean_key in valid_keys:
                    clean_key = valid_keys[clean_key]
                else:
                    invalid_keys.append(key)
                    break
                    
                if clean_key in path_keys:
                    value = setup_absolute_paths(path_to_file,value)
                configdict[valid_keys[clean_key]] = value
                overwriting += 1

    if len(invalid_keys)==1:
        sys.stderr.write(cyan(f'Error: invalid key in config file.\n') + f'\t- {invalid_keys[0]}\n')
        sys.exit(-1)
    elif len(invalid_keys) >1:
        keys = ""
        for i in invalid_keys:
            keys += f"\t- {i}\n"
        sys.stderr.write(cyan(f'Error: invalid keys in config file.\n') + f'{keys}')
        sys.exit(-1)
    print(green(f"Adding {overwriting} arguments to internal config."))

def setup_config_dict(cwd,config_arg):
    config = get_defaults()
    config["cwd"] = cwd
    if config_arg:
        configfile = check_configfile(cwd,config_arg)
        
        parse_yaml_file(configfile,config)

    else:
        config["input_path"] = cwd
    return config

def misc_args_to_config(verbose,threads, civet_mode, config):
    misc.add_arg_to_config("verbose",verbose,config)
    misc.add_arg_to_config("threads",threads,config)
    misc.add_arg_to_config("civet_mode", civet_mode, config)


def set_up_verbosity(config):
    if config["verbose"]:
        config["quiet"] = False
        config["log_api"] = ""
        config["log_string"] = ""
    else:
        config["quiet"] = True
        logger = custom_logger.Logger()
        config["log_api"] = logger.log_handler

        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "