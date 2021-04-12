#!/usr/bin/env python3
from civet.utils import log_colours as colour
from civet.utils import misc
import sys
import os
import csv
from Bio import SeqIO


"""
i_group parsing
- Either csv or id string
- Is the id string comma separated, unique the ids
- Does the csv file exist, is it a .csv file, is it readable by the csv module
- is the fasta file in fasta format
- are maxambig and minlen the correct format in the right range
"""
def ids_qc(ids):
    if not type(ids) is list:
        ids = ids.split(",")
        
    ids = list(set(ids))
    print(colour.green("Unique IDs input: ") +  f"{len(ids)}")
    return ids

def csv_qc(input_csv,input_column):

    ending = input_csv.split(".")[-1]

    if ending in ["yaml","yml","json"]:
        sys.stderr.write(colour.cyan(f"Error: -i,--input-csv accepts a csv file. As of civet 3.0 please use -c/ --config to input a config file or -ids/ --id-string to input a comma-separated string of IDs.\n"))
        sys.exit(-1)
    elif ending in ["xls","xlsx"]:
        sys.stderr.write(colour.cyan(f"Error: it looks like you've provided an excel file as input.\n-i,--input-csv accepts a csv file\n"))
        sys.exit(-1)
    elif ending == "csv":
        pass
    else:
        sys.stderr.write(colour.cyan(f"Error: -i,--input-csv accepts a csv file. As of civet 3.0 please use -c/ --config to input a config file or -ids/ --id-string to input a comma-separated string of IDs.\n"))
        sys.exit(-1)
    
    if os.path.isfile(input_csv):
        try:
            c = 0
            with open(input_csv,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    c +=1
            print(colour.green(f"{c} rows in input csv file."))
        except:
            sys.stderr.write(colour.cyan(f"Unable to read csv file, please check your input csv is in the correct format with a header.\n"))
            sys.exit(-1)
    else:
        sys.stderr.write(colour.cyan(f"Cannot find csv file: ")+f"{input_csv}\n")
        sys.exit(-1)
    print(colour.green(f"Input csv file:") + f" {input_csv}")

    with open(input_csv,"r") as f:
        reader = csv.DictReader(f)
        if input_column in reader.fieldnames:
            pass
        else:
            encode = False
            for i in reader.fieldnames:
                if "ufeff" in i:
                    sys.stderr.write(colour.cyan(f"Error: it appears your csv file may have been edited in Excel and now contains hidden characters.\n") + "Please remove said characters in a text editor and try again.")
                    sys.exit(-1)
            else:
                sys.stderr.write(colour.cyan(f"Error: {input_column} column not found in input csv file.\n"))
                sys.exit(-1)

def input_query_parsing(input_csv,input_column,ids,config):

    misc.add_arg_to_config("ids",ids,config)

    misc.add_file_to_config("input_csv",input_csv,config)
    misc.add_arg_to_config("input_column",input_column,config)

    if "ids" in config and "input_csv" in config:
        sys.stderr.write(colour.cyan(f"Error: it looks like you've provide a csv file and an ID string, please provide one or the other.\n"))
        sys.exit(-1)
    elif "ids" in config:
        config["ids"] = ids_qc(config["ids"])
        
    elif "input_csv" in config:
        csv_qc(config["input_csv"],config["input_column"])
    


def input_fasta_check(input_file):

    ending = input_file.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(colour.cyan(f"Please input sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)
    else:
        print(colour.green(f"Input fasta file:") + f" {input_file}")
    c = 0
    for record in SeqIO.parse(input_file, "fasta"):
        c += 1
        
    if c == 0:
        sys.stderr.write(colour.cyan(f"Error: no records found in fasta file, please check your input fasta is in the correct format.\n"))
        sys.exit(-1)
    else:
        print(colour.green(f"{c} records in input fasta file."))


def fasta_qc_level(maxambig,minlen):
    if not type(minlen) is int:
        try:
            minlen = int(minlen)
        except:
            sys.stderr.write(colour.cyan(f"Error: --min-length must be a positive integer.\n"))
            sys.exit(-1)
    if not minlen > 0:
        sys.stderr.write(colour.cyan(f"Error: --min-length must be a positive integer.\n"))
        sys.exit(-1)
    
    if not type(maxambig) is float:
        try:
            maxambig = float(maxambig)
        except:
            sys.stderr.write(colour.cyan(f"Error: --max-ambiguity must be between 0 and 1.\n"))
            sys.exit(-1)
    if not (maxambig <=1 and maxambig >=0):
        sys.stderr.write(colour.cyan(f"Error: --max-ambiguity must be between 0 and 1.\n"))
        sys.exit(-1)

def fasta_ids_list(fasta):
    ids = []
    for record in SeqIO.parse(fasta,"fasta"):
        ids.append(record.id)
    ids = ids_qc(ids)

    return ids

def input_fasta_parsing(input_fasta,maxambig,minlen,config):
    
    misc.add_file_to_config("fasta",input_fasta,config)
    misc.add_arg_to_config("max_ambiguity",maxambig,config)
    misc.add_arg_to_config("min_length",minlen,config)

    if "fasta" in config:
        
        input_fasta_check(config["fasta"])
        fasta_qc_level(config["max_ambiguity"],config["min_length"])

        if not "input_csv" in config and not "ids" in config:
            config["ids"] = fasta_ids_list(config["fasta"])
