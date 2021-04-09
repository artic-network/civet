from datetime import date
import os
import sys
import yaml
from civet.utils import log_colours as colour

def get_defaults():
    today = date.today()
    default_dict = {
                    "date": today,# date investigation was opened
                    "output_prefix":"civet",
                    "authors": "", # List of authors, affiliations and contact details

                    # Initialising variables
                    "num_seqs":0,
                    "datadir":"",
                    
                    # Search defined by metadata column
                    "from_metadata":False,

                    # Columns to match
                    "input_column":"name",
                    "database_column":"central_sample_id",
                    "database_date_column":"sample_date",
                    "input_date_column":"sample_date",

                    # Search parameters
                    "distance":2,
                    "up_distance":False,
                    "down_distance":False,
                    "collapse_threshold":1,
                    
                    # reference coordinates to pad until and from
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674,

                    # QC standards for input fasta file
                    "max_ambiguity":0.5,
                    "min_length":10000,

                    # misc defaults
                    "threads":1,
                    "force":True,
                    "no_temp":False,
                    }
    return default_dict

def check_configfile(cwd,config_arg):
    configfile = os.path.join(cwd,config_arg)

    ending = configfile.split(".")[-1]

    if ending not in ["yaml","yml"]:
        sys.stderr.write(colour.cyan(f'Error: config file {configfile} must be in yaml format.\n'))
        sys.exit(-1)
    
    elif not os.path.isfile(configfile):
        sys.stderr.write(colour.cyan(f'Error: cannot find config file at {configfile}\n'))
        sys.exit(-1)
    else:
        print(colour.green(f"Input config file:") + f" {configfile}")
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
                # ogroup args

                # misc group args
                "t":"threads",
                "threads":"threads"
    }
    for i in config:
        arguments[i] = i
    return arguments

def load_yaml(f):
    try:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    except:
        sys.stderr.write(colour.cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
        sys.exit(-1)
    return input_config

def setup_file_paths(config):
    input_files = ["input_csv","fasta","background_csv","background_fasta"]

    for key in config:
        if key in input_files:
            config[key] = os.path.join(config["input_path"],config[key])


def parse_yaml_file(configfile,configdict):
    overwriting = 0

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
                    configdict[valid_keys[clean_key]] = value
                    overwriting += 1
                else:
                    invalid_keys.append(key)

    if len(invalid_keys)==1:
        sys.stderr.write(colour.cyan(f'Error: invalid key in config file.\n') + f'\t- {invalid_keys[0]}\n')
        sys.exit(-1)
    elif len(invalid_keys) >1:
        keys = ""
        for i in invalid_keys:
            keys += f"\t- {i}\n"
        sys.stderr.write(colour.cyan(f'Error: invalid keys in config file.\n') + f'{keys}')
        sys.exit(-1)
    print(colour.green(f"Adding {overwriting} arguments to internal config."))

def setup_config_dict(cwd,config_arg):
    config = get_defaults()
    config["cwd"] = cwd
    if config_arg:
        configfile = check_configfile(cwd,config_arg)
        
        parse_yaml_file(configfile,config)

        setup_file_paths(config)
    else:
        config["input_path"] = cwd
    return config
