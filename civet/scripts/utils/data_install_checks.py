#!/usr/bin/env python3
import pkg_resources
import log_colours as colour
import sys
import os


def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('civet', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(colour.cyan(f'Error: Missing package data.')+f'\n\t- {filename}\nPlease install the latest civet version with `civet --update`.\n')
        sys.exit(-1)
    
def check_install():
    resources = [
        {"key":"reference_fasta",
        "directory":"data",
        "filename":"reference.fasta"},
        {"key":"outgroup_fasta",
        "directory":"data",
        "filename":"outgroup.fasta"}
    ]
    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)

# config={}
# check_install()
