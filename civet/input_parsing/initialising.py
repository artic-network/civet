from datetime import date
import os
import sys
import yaml
from civet.utils.log_colours import green,cyan

import civet.utils.custom_logger as custom_logger
import civet.utils.log_handler_handle as lh
from civet.utils import misc
from civet.utils.config import *

def get_defaults():
    today = date.today()
    default_dict = {
                    
                    #input options
                    KEY_FROM_METADATA:False,
                    KEY_MAX_QUERIES:5000,

                    #input cols
                    KEY_INPUT_ID_COLUMN:"name",
                    KEY_INPUT_DATE_COLUMN:False, #default is sample_date if present in the input csv
                    KEY_INPUT_DISPLAY_COLUMN: False,

                    # input seq options 
                    KEY_NUM_SEQS:0,
                    KEY_MAX_AMBIGUITY:0.5,
                    KEY_MIN_LENGTH:20000,
                    KEY_TRIM_START:265,   # where to pad to using datafunk
                    KEY_TRIM_END:29674,

                    # background variables
                    
                    CIVET_DATADIR: os.getenv('CIVET_DATADIR'),
                    KEY_BACKGROUND_METADATA:False,
                    KEY_BACKGROUND_SEQUENCES:False,
                    KEY_BACKGROUND_SNPS:False,
                    KEY_BACKGROUND_TREE:False,

                    # background columns
                    
                    KEY_BACKGROUND_ID_COLUMN:"sequence_name",
                    KEY_SEQUENCE_ID_COLUMN:False,
                    KEY_BACKGROUND_DATE_COLUMN:False, #default is sample_date if present in the background csv, then date_column if present
                    KEY_BACKGROUND_LOCATION_COLUMN:False, #default is country if present in the background csv

                    # Output options
                    "output_prefix":"civet",
                    "output_data":False,
                    "datestamp":False,
                    "no_temp":False,
                    "overwrite":False,

                    # catchment config
                    "catchment_background_size":100,
                    "snp_distance":False,
                    "snp_distance_up":2,
                    "snp_distance_down":2,
                    "snp_distance_side":2,
                    "collapse_threshold":1,

                    "downsample":["mode=random"],
                    "mode":False,
                    "factor":False,
                    "downsample_field":False,
                    "downsample_column":False,                    

                    #Report options
                    "colour_theme": "#7178bc",
                    "colour_map": "#86b0a6,#565b7b,#e9c46a,#e1998a,#d19c2c,#264653,#4888a3,#c77939,#b56576,#eb7e83,#F46D43,#9e2a2b,#84a98c",
                    KEY_REPORT_CONTENT:["1,2,3,4,5"],
                    KEY_REPORT_PRESET:False,
                    "report_title":f"civet report",
                    KEY_ANONYMISE:False,
                    
                    #Table options
                    "query_table_content":False,
                    KEY_MUTATIONS:[],

                    #Tree options
                    "tree_annotations": "lineage,country",
                    "max_tree_size":500,

                    #Timeline options
                    "timeline_dates":False, #default is date_column or background_date_column if they are present
                    "timeline_group_column":False, #default is input_display_column

                    #Background map options
                    "background_map_column":False, #default is location_column
                    "background_map_location":False, #default is all valid locations in background metadata
                    "background_map_file": False, #Same as above
                    "centroid_file":False, #will select appropriate centroid file for provided jsons as default

                    #Query map options
                    KEY_QUERY_MAP_FILE:False, #if UK, will find the geojson containing aggregated_adm2s, if global will find adm0 or adm1 geojson
                    KEY_LATITUDE_COLUMN:"latitude",
                    KEY_LONGITUDE_COLUMN:"longitude",
                    
                    # "map_column":False,
                    KEY_BACKGROUND_MAP_DATE_RANGE:False, #default is no restriction, will do background from 2019-11-01 to today

                    # misc defaults
                    KEY_CIVET_MODE: os.getenv('CIVETMODE'),
                    KEY_THREADS:1,
                    KEY_VERBOSE:False,
                    KEY_DATE: today,# date investigation was opened
                    KEY_AUTHORS: "" # List of authors, affiliations and contact details
                    
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
                "ids":"id_string",
                "id_string":"id_string",
                "i":"input_metadata",
                "input_metadata":"input_metadata",
                "f":"input_sequences",
                "input_sequences":"input_sequences",
                "fm":"from_metadata",
                "from_metadata":"from_metadata",
                "max_queries":"max_queries",
                "mq":"max_queries",

                # ic group args
                "icol":"input_id_column",
                "input_id_column":"input_id_column",
                "input_date_column":"input_date_column",
                "idate":"input_date_column",
                "idisp":"input_display_column",
                "input_display_column":"input_display_column",

                # is group args
                "n":"max_ambiguity",
                "max_ambiguity":"max_ambiguity",
                "l":"min_length",
                "min_length":"min_length",
                "ts":"trim_start",
                "te":"trim_end",
                "trim_start":"trim_start",
                "trim_end":"trim_end",

                # dgroup args
                "datadir":"datadir",
                "DATADIR":"datadir",
                "d":"datadir",
                "background_metadata":"background_metadata",
                "bm":"background_metadata",
                "background_sequences":"background_sequences",
                "bseq":"background_sequences",
                "background_snps":"background_snps",
                "bsnp":"background_snps",
                "background_tree":"background_tree",
                "bt":"background_tree",

                #bc group args
                "bicol":"background_id_column",
                "background_id_column":"background_id_column",
                "sicol":"sequence_id_column",
                "sequence_id_column":"sequence_id_column",
                "background_date_column":"background_date_column",
                "bdate":"background_date_column",
                "background_location_column":"background_location_column",
                "bloc":"background_location_column",

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

                # cgroup args
                
                "cs":"catchment_background_size",
                "catchment_background_size":"catchment_background_size",
                "downsample":"downsample",
                "ds":"downsample",
                "snp_distance":"snp_distance",
                "snpd":"snp_distance",
                "snp_distance_up":"snp_distance_up",
                "snp_distance_down":"snp_distance_down",
                "snp_distance_side":"snp_distance_side",

                #rgroup args
                "report_content":"report_content",
                "rc":"report_content",
                "report_preset":"report_preset",
                "rp":"report_preset",
                "report_title":"report_title",
                "rt":"report_title",
                "colour_theme":"colour_theme",
                "ct":"colour_theme",
                "colour_map":"colour_map",
                "cmap":"colour_map",
                
                "anonymise":"anonymise",
                "anonymize":"anonymise",

                #tb group args
                "query_table_content":"query_table_content",
                "catchment_table_content":"catchment_table_content",
                "mutations":"mutations",

                #t group args
                "tree_annotations":"tree_annotations",
                "ta":"tree_annotations",
                "max_tree_size":"max_tree_size",
                "mt":"max_tree_size",

                #tl group args
                "td":"timeline_dates",
                "timeline_dates":"timeline_dates",
                "timeline_group_column":"timeline_group_column",
                "tgc":"timeline_group_column",
                
                #bm group args
                "background_map_date_range":"background_map_date_range",
                "bmdr":"background_map_date_range",
                "background_map_column":"background_map_column",
                "bmcol":"background_map_column",
                "background_map_file":"background_map_file",
                "background_map_location":"background_map_location",
                "bmloc":"background_map_location",
                "bmfile":"background_map_file",
                "centroid_file":"centroid_file",

                #qm group args
                "lat":"latitude_column",
                "latitude_column":"latitude_column",
                "long":"longitude_column",
                "longitude_column":"longitude_column",
                "query_map_file":"query_map_file",
                "qmfile":"query_map_file",
                
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
    return [KEY_INPUT_METADATA,KEY_INPUT_SEQUENCES,KEY_BACKGROUND_METADATA,KEY_BACKGROUND_SEQUENCES,KEY_BACKGROUND_TREE,KEY_BACKGROUND_SNPS,CIVET_DATADIR,KEY_OUTDIR,KEY_TEMPDIR]

def setup_absolute_paths(path_to_file,value):
    return os.path.join(path_to_file,value)


def parse_yaml_file(configfile,configdict):
    overwriting = 0
    path_keys = return_path_keys()

    path_to_file = os.path.abspath(os.path.dirname(configfile))
    configdict[KEY_INPUT_PATH] = path_to_file
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
    config[KEY_CWD] = cwd
    if config_arg:
        configfile = check_configfile(cwd,config_arg)
        
        parse_yaml_file(configfile,config)

    else:
        config[KEY_INPUT_PATH] = cwd
    return config

def misc_args_to_config(verbose,threads, civet_mode, config):
    misc.add_arg_to_config(KEY_VERBOSE,verbose,config)
    misc.add_arg_to_config(KEY_THREADS,threads,config)
    misc.add_arg_to_config(KEY_CIVET_MODE, civet_mode, config)


def set_up_verbosity(config):
    if config[KEY_VERBOSE]:
        config[KEY_QUIET] = False
        config[KEY_LOG_API] = ""
        config[KEY_LOG_STRING] = ""
    else:
        config[KEY_QUIET] = True
        logger = custom_logger.Logger()
        config[KEY_LOG_API] = logger.log_handler

        lh_path = os.path.realpath(lh.__file__)
        config[KEY_LOG_STRING] = f"--quiet --log-handler-script {lh_path} "