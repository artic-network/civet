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
        "directory":"data",
        "filename":"report_template.mako"},
        {"key":"query_chunk",
        "directory":"data/report_chunks",
        "filename":"query_table_chunk.txt"},
        {"key":"timeline_chunk",
        "directory":"data/report_chunks",
        "filename":"timeline_chunk.txt"}
    ]
    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)

# config={}
# check_install()
