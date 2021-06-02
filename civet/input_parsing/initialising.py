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
                    "data_column":"sequence_name",
                    "fasta_column":False,

                    # Search parameters
                    "distance":2,
                    "up_distance":False,
                    "down_distance":False,
                    "collapse_threshold":1,
                    
                    # reference coordinates to pad until and from
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674,

                    # catchment config
                    "catchment_size":1000,

                    # QC standards for input fasta file
                    "max_ambiguity":0.5,
                    "min_length":20000,

                    #Report options
                    "alt_seq_name":"sequence_name"
                    "anonymise":False,
                    "timeline_dates":False,
                    "found_table_cols":False,
                    "provided_table_cols":False,
                    "data_date_column":"sample_date",
                    "input_date_column":"sample_date",

                    #Map options
                    "map_queries":False,
                    "latitude_column":"latitude",
                    "longitude_column":"longitude",
                    "map_background":False,
                    "background_map_column":"adm1",
                    "background_map_date_window":7,
                    "background_map_date_start":False,
                    "background_map_date_end":False,
                    "background_map_date_column":"sample_date",

                    # misc defaults
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
                "dcol":"data_column",
                "data_column":"data_column",
                "fcol":"fasta_column",
                "fasta_column":"fasta_column",

                # ogroup args
                "o":"outdir",
                "outdir":"outdir",
                "p":"output_prefix",
                "output_prefix":"output_prefix",
                "datestamp":"datestamp",
                "ds":"datestamp",
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

                #rgroup args
                "alt":"alt_seq_name",
                "alt_seq_name":"alt_seq_name"
                "anon":"anonymise",
                "anonymise":"anonymise",
                "anonymize":"anonymise",
                "td":"timeline_dates",
                "timeline_dates":"timeline_dates",
                "ftcol":"found_table_cols",
                "found_table_cols":"found_table_cols",
                "ptcol":"provided_table_cols".
                "provided_table_cols":"provided_table_cols",
                "data_date_column":"data_date_column",
                "ddc":"data_date_column",
                "input_date_column":"input_date_column",
                "idatc":"input_date_column",

                #mgroup args
                "mq":"map_queries",
                "map_queries":"map_queries"
                "lat":"latitude_column",
                "latitude_column":"latitude_column",
                "long":"longitude_column",
                "longitude_column":"longitude_column"
                "mbg":"map_background",
                "map_background":"map_background",
                "bgcol":"background_map_column",
                "background_map_column":"background_map_column",
                "dw":"background_map_date_window",
                "background_map_date_window":"background_map_date_window",
                "ds":"background_map_date_start",
                "background_map_date_start":"background_map_date_start",
                "de":"background_map_date_end",
                "background_map_date_end":"background_map_date_end",
                "bgdate":"background_map_date_column",
                "background_map_date_column":"background_map_date_column",

                # misc group args
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

def misc_args_to_config(verbose,threads,config):
    misc.add_arg_to_config("verbose",verbose,config)
    misc.add_arg_to_config("threads",threads,config)

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