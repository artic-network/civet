#!/usr/bin/env python3
import pkg_resources
from civet.utils.log_colours import green,cyan
import sys
import os


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
    resources = [
        {"key":"reference_fasta",
        "directory":"data",
        "filename":"reference.fasta"},
        {"key":"outgroup_fasta",
        "directory":"data",
        "filename":"outgroup.fasta"},
        {"key":"report_template",
        "directory":"data/report_modules",
        "filename":"report_template.mako"},
        {"key":"mako_query_table",
        "directory":"data/report_modules",
        "filename":"query_table.txt"},
        {"key":"mako_timeline",
        "directory":"data/report_modules",
        "filename":"timeline.txt"},
        {"key":"html_colours",
        "directory":"data",
        "filename":"html_colours.csv"},
        {"key":"global_accepted_values",
        "directory":"data/map_data",
        "filename":"global_accepted_values.tsv"},
        {"key":"uk_acceptable_values",
        "directory":"data/map_data",
        "filename":"UK_acceptable_values.tsv"}
    ]
    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)

# config={}
# check_install()
