#!/usr/bin/env python3
import subprocess
import os
import sys
import log_colours as colour
import importlib

def which(dependency):
    try:
        subprocess.check_output(["which", dependency])
        return True
    except subprocess.CalledProcessError:
        return False

def check_module(module, missing):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)

def check_this_dependency(dependency,missing):
    check = which(dependency)

    if not check:
        missing.append(dependency)

def check_dependencies():

    missing = []

    dependency_list = ["gofasta","minimap2","snakemake","iqtree","jclusterfunk"]
    module_list = ["baltic","mako","Bio"]

    for dependency in dependency_list:
        check_this_dependency(dependency, missing)

    for module in module_list:
        check_module(module, missing)

    if missing:
        if len(missing)==1:
            sys.stderr.write(colour.cyan(f'Error: Missing dependency `{missing[0]}`.')+'\nPlease update your civet environment.\n')
            sys.exit(-1)
        else:
            dependencies = ""
            for i in missing:
                dependencies+=f"\t- {i}\n"

            sys.stderr.write(colour.cyan(f'Error: Missing dependencies.')+f'\n{dependencies}Please update your civet environment.\n')
            sys.exit(-1)
    else:
        print(colour.green("All dependencies satisfied."))


# check_dependencies()
