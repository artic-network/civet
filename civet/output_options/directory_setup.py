#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
from civet.utils import misc

import sys
import os
import glob

from datetime import datetime 
from datetime import date
import tempfile
import sys

"""
Desired behaviour

Default outdir -> civet-2021-XX-YY,
check if that exists, append a number if it already exists -> civet-2021-XX-YY-2
--overwrite flag to overwrite the civet directory

If -o/--outdir then output that as the outdir instead, without a datestamp
If -p/--output-prefix then output directory as: prefix-2021-XX-YY

prefix-2021-XX-YY (or civet-2021-XX-YY) will also be the name of the final civet report.

"""

def datestamped_outdir(config):
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    if not "outdir" in config:
        expanded_path = os.path.expanduser(config["cwd"])
        outdir = os.path.join(expanded_path, f"{config['output_prefix']}_{d}")
    else:
        if config["datestamp"]:
            outdir = f"{config['outdir']}_{d}"
        else:
            outdir = config['outdir']

    if not config["overwrite"]:
        counter = 1
        while os.path.exists(outdir):
            if outdir.split("_")[-1].isdigit():
                outdir = "_".join(outdir.split("_")[:-1])
                outdir = f"{outdir}_{counter}"
            else:
                outdir = f"{outdir}_{counter}"
            counter +=1

    return outdir,d

def clear_old_files(config):
    if config["overwrite"] and os.path.exists(config["outdir"]):
        print(green("Overwriting previous output in ") + config["outdir"] + ".")
        old_files = glob.glob(f'{config["outdir"]}/*.*', recursive=True)

        for f in old_files:
            try:
                os.remove(f)
                print(green("Removed: "),f)
            except:
                print(cyan("Can't remove "),f)

def output_report_filenames(d,config):

    reports = config["reports"]
    report_names = []
    c = 0
    if len(reports)>1:
        for i in reports:
            c +=1
            if config["datestamp"]:
                output_report = f"{config['output_prefix']}_{d}_{c}.html"
            else:
                output_report = f"{config['output_prefix']}_{c}.html"
            report_names.append(output_report)
    else:
        if config["datestamp"]:
            output_report = f"{config['output_prefix']}_{d}.html"
        else:
            output_report = f"{config['output_prefix']}.html"
        report_names.append(output_report)
    return report_names

def set_up_tempdir(config):
    if config["no_temp"]:
        tempdir = config["outdir"]
        config["tempdir"] = tempdir
    elif "tempdir" in config:
        to_be_dir = config["tempdir"]
        try:
            if not os.path.exists(to_be_dir):
                os.mkdir(to_be_dir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {to_be_dir}.\n'))
            sys.exit(-1)
        tempdir = tempfile.mkdtemp(dir=to_be_dir)
        config["tempdir"] = tempdir
    else:
        tempdir = tempfile.mkdtemp()
        config["tempdir"] = tempdir
        try:
            if not os.path.exists(tempdir):
                os.mkdir(tempdir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {tempdir}.\n'))
            sys.exit(-1)
        
        try:
            with open(os.path.join(tempdir, "test.txt"),"w") as fw:
                fw.write("Test")
        except:
            sys.stderr.write(cyan(f'Error: cannot write to temp directory {tempdir}.\n'))
            sys.exit(-1)

def set_up_data_outdir(config):
    if config["output_data"]:
        config["data_outdir"] = config["outdir"]
    else:
        config["data_outdir"] = config["tempdir"]

def output_group_parsing(outdir,output_prefix,overwrite,datestamp,output_data,tempdir,no_temp,config):
    
    misc.add_path_to_config("outdir",outdir,config)
    misc.add_arg_to_config("output_prefix",output_prefix,config)
    misc.add_arg_to_config("overwrite",overwrite,config)
    misc.add_arg_to_config("datestamp",datestamp,config)
    misc.add_arg_to_config("output_data",output_data,config)
    misc.add_path_to_config("tempdir",tempdir,config)
    misc.add_arg_to_config("no_temp",no_temp,config)

    config["outdir"],d = datestamped_outdir(config)
    
    set_up_tempdir(config)
    set_up_data_outdir(config)

    config["output_reports"] = output_report_filenames(d,config)

    clear_old_files(config)

    if not os.path.exists(config["outdir"]):
        os.mkdir(config["outdir"])




