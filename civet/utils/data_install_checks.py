#!/usr/bin/env python3
import pkg_resources
from civet.utils.log_colours import green,cyan
import sys
import os
from civet.utils.config import *


def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('civet', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(cyan(f'Error: Missing package data.')+f'\n\t- {filename}\nPlease install the latest civet version with `civet --update`.\n')
        sys.exit(-1)

def get_snakefile(thisdir):
    snakefile = os.path.join(thisdir, 'scripts','civet.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n'))
        sys.exit(-1)
    return snakefile

def check_install(config):
    for resource in resources:
        package_data_check(resource[RESOURCE_KEY_FILENAME],resource[RESOURCE_KEY_DIRECTORY],resource[RESOURCE_KEY],config)

# config={}
# check_install()
