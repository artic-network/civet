from datetime import date
import os
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
                    "input_id_column":"name",
                    "database_id_column":"central_sample_id",
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



def parse_yaml_file(configfile,configdict):
    overwriting = 0

    with open(configfile,"r") as f:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
        for key in input_config:
            snakecase_key = key.replace("-","_").lstrip("-")
            configdict[snakecase_key] = input_config[key]
            overwriting += 1
    print(colour.green(f"Adding {overwriting} arguments to internal config."))


def setup_config_dict(cwd,config_arg):
    config = get_defaults()

    if config_arg:
        configfile = check_configfile(cwd,config_arg)
        
        parse_yaml_file(configfile,config)

    return config



