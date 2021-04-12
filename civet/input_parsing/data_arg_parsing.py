#!/usr/bin/env python3
from civet.utils import log_colours as colour
from civet.utils import misc
import sys
import os
import csv
from Bio import SeqIO

def check_dir_for_file(file_list, extension,argument,key,config):
    match = []
    for fn in file_list:
        if fn.endswith(extension):
            match.append(fn)
    if len(match) > 1:
        sys.stderr.write(colour.cyan(f"Error: More than 1 {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
        sys.exit(-1)
    elif len(match) == 0:
        sys.stderr.write(colour.cyan(f"Error: No {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
        sys.exit(-1)
    elif len(match) == 1:
        config[key] = match[0]

def background_data_load(config):
    file_list = []
    for r,d,f in os.walk(config["datadir"]):
        for fn in f:
            file_list.append(os.path.join(r,fn))
    if not config["background_csv"]:
        check_dir_for_file(file_list, "csv","-bc/--background-csv","background_csv",config)
    if not config["background_fasta"]:
        check_dir_for_file(file_list, "fasta","-bf/--background-fasta","background_fasta",config)

def check_datadir(config):
    datadir = config["datadir"]
    if not os.path.exists(datadir):
        if not config["background_csv"] or not config["background_fasta"]:
            print(colour.cyan(f"Error: data directory not found: {datadir}.\n")+"Please check the datadir path (-d/ --datadir) provided exists or supply files with -bc/ --background-csv and -bf/ --background-fasta.\n")
            sys.exit(-1)

    background_data_load(config)


def check_background_csv(background_csv,data_column):
    
    ending = background_csv.split(".")[-1]

    if ending in ["xls","xlsx"]:
        sys.stderr.write(colour.cyan(f"Error: it looks like you've provided an excel file as input.\n-bc,--background-csv accepts a csv file.\n"))
        sys.exit(-1)
    elif ending == "csv":
        pass
    else:
        sys.stderr.write(colour.cyan(f"Error: please provide a csv file for -bc/ --background-csv.\n"))
        sys.exit(-1)
    
    if os.path.isfile(background_csv):
        try:
            c = 0
            with open(background_csv,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    c +=1
            print(colour.green(f"{c} rows in background csv file."))
        except:
            sys.stderr.write(colour.cyan(f"Unable to read background csv file, please check your background csv is in the correct format with a header.\n"))
            sys.exit(-1)
    else:
        sys.stderr.write(colour.cyan(f"Cannot find csv file: ")+f"{background_csv}\n")
        sys.exit(-1)
    print(colour.green(f"Background csv file:") + f" {background_csv}")

    with open(background_csv,"r") as f:
        reader = csv.DictReader(f)
        if data_column in reader.fieldnames:
            pass
        else:
            encode = False
            for i in reader.fieldnames:
                if "ufeff" in i:
                    sys.stderr.write(colour.cyan(f"Error: it appears your csv file may have been edited in Excel and now contains hidden characters.\n") + "Please remove said characters in a text editor and try again.")
                    sys.exit(-1)
            else:
                sys.stderr.write(colour.cyan(f"Error: {data_column} column not found in background csv file. Please specifiy which column to match with -dcol/ --data-column.\n"))
                sys.exit(-1)
    return c

def check_background_fasta(background_fasta):
    ending = background_fasta.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(colour.cyan(f"Please input background sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)
    else:
        print(colour.green(f"Background fasta file:") + f" {background_fasta}")
    c = 0
    for record in SeqIO.parse(background_fasta, "fasta"):
        c += 1
        
    if c == 0:
        sys.stderr.write(colour.cyan(f"Error: no records found in fasta file, please check your background fasta is in the correct format.\n"))
        sys.exit(-1)
    else:
        print(colour.green(f"{c} records in background fasta file."))
    
    return c

def data_group_parsing(datadir,background_csv,background_fasta,data_column,config):
    
    if datadir:
        datadir = os.path.abspath(datadir)

    misc.add_arg_to_config("datadir",datadir,config)
    misc.add_file_to_config("background_csv",background_csv,config)
    misc.add_file_to_config("background_fasta",background_fasta,config)

    if not config["datadir"]:
        if not config["background_csv"] and not config["background_fasta"]:
            sys.stderr.write(colour.cyan(f"Error: insufficient background data supplied, please provide a background csv and fasta file.\n"))
            sys.exit(-1)
        
    if config["datadir"]:
        check_datadir(config)
    
    misc.add_arg_to_config("data_column",data_column,config)
    
    background_data_load(config)
    
    csv_record_count = check_background_csv(config["background_csv"],config["data_column"])
    fasta_record_count = check_background_fasta(config["background_fasta"])
    
    if csv_record_count != fasta_record_count:
        sys.stderr.write(colour.cyan(f"Error: different number of background csv and background fasta records.\n")+"Please provide a fasta record for each row in the background metadata file.\n")
        sys.exit(-1)
