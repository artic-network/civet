#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import sys
import os
import csv
from Bio import SeqIO

from datetime import datetime 
from datetime import date
from civet.utils.config import *


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
    print(green("Unique IDs input: ") +  f"{len(ids)}")
    return ids

def check_for_protected_col_names(header):
    for field in header:
        if field in PROTECTED_COL_NAMES:
            sys.stderr.write(cyan(f"Error: `{field}` is a protected column name used internally in civet, please rename this column.\n"))
            sys.exit(-1)

def csv_qc(input_metadata,input_id_column):

    ending = input_metadata.split(".")[-1]

    if ending in ["yaml","yml","json"]:
        sys.stderr.write(cyan(f"Error: -i,--input-metadata accepts a csv file. As of civet 3.0 please use -c/ --config to input a config file or -ids/ --id-string to input a comma-separated string of IDs.\n"))
        sys.exit(-1)
    elif ending in ["xls","xlsx"]:
        sys.stderr.write(cyan(f"Error: it looks like you've provided an excel file as input.\n-i,--input-metadata accepts a csv file\n"))
        sys.exit(-1)
    elif ending == "csv":
        pass
    else:
        sys.stderr.write(cyan(f"Error: -i,--input-metadata accepts a csv file. As of civet 3.0 please use -c/ --config to input a config file or -ids/ --id-string to input a comma-separated string of IDs.\n"))
        sys.exit(-1)
    
    if os.path.isfile(input_metadata):
        try:
            c = 0
            with open(input_metadata,"r") as f:
                reader = misc.read_csv_or_tsv(input_metadata,f)
                for row in reader:
                    c +=1
            # print(green(f"{c} rows in input csv file."))
        except:
            sys.stderr.write(cyan(f"Unable to read csv file, please check your input csv is in the correct format with a header.\n"))
            sys.exit(-1)
    else:
        sys.stderr.write(cyan(f"Cannot find csv file: ")+f"{input_metadata}\n")
        sys.exit(-1)
    print(green(f"Input csv file:") + f" {input_metadata}.")

    input_ids = []
    with open(input_metadata,"r") as f:
        reader = misc.read_csv_or_tsv(input_metadata,f)
        if input_id_column in reader.fieldnames:
            for row in reader:
                input_ids.append(row[input_id_column])
        else:
            encode = False
            for i in reader.fieldnames:
                if "ufeff" in i:
                    sys.stderr.write(cyan(f"Error: it appears your csv file may have been edited in Excel and now contains hidden characters.\n") + "Please remove said characters in a text editor and try again.")
                    sys.exit(-1)
            else:
                sys.stderr.write(cyan(f"Error: {input_id_column} column not found in input csv file.\n"))
                sys.exit(-1)
    return input_ids

def input_query_parsing(input_metadata,input_id_column,ids,from_metadata,config):

    misc.add_arg_to_config(KEY_IDS,ids,config)

    misc.add_file_to_config(KEY_INPUT_METADATA,input_metadata,config)
    misc.add_arg_to_config(KEY_INPUT_ID_COLUMN,input_id_column,config)

    misc.add_arg_to_config(KEY_FROM_METADATA,from_metadata,config)

    if KEY_IDS in config and KEY_INPUT_METADATA in config:
        sys.stderr.write(cyan(f"Error: it looks like you've provided a csv file and an ID string, please provide one or the other.\n"))
        sys.exit(-1)
    elif KEY_IDS in config:
        config[KEY_IDS] = ids_qc(config[KEY_IDS])
        
    elif KEY_INPUT_METADATA in config:
        config[KEY_IDS] = csv_qc(config[KEY_INPUT_METADATA],config[KEY_INPUT_ID_COLUMN])

    if config[KEY_FROM_METADATA] and KEY_IDS in config:
        sys.stderr.write(cyan('Error: civet accepts either -fm/--from-metadata (which generates a query from the background data) or an input query (supplied via `-ids/--id-string`, `-i/--input-metadata` or `-f/--input-sequences`).\n'))
        sys.exit(-1)
    
def input_fasta_check(input_file):

    ending = input_file.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(cyan(f"Please input sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)
    else:
        print(green(f"Input fasta file:") + f" {input_file}")
    c = 0
    for record in SeqIO.parse(input_file, "fasta"):
        c += 1
        
    if c == 0:
        sys.stderr.write(cyan(f"Error: no records found in fasta file, please check your input fasta is in the correct format.\n"))
        sys.exit(-1)
    # else:
    #     print(green(f"{c} records in input fasta file."))


def fasta_qc_level(maxambig,minlen):
    if not type(minlen) is int:
        try:
            minlen = int(minlen)
        except:
            sys.stderr.write(cyan(f"Error: --min-length must be a positive integer.\n"))
            sys.exit(-1)
    if not minlen > 0:
        sys.stderr.write(cyan(f"Error: --min-length must be a positive integer.\n"))
        sys.exit(-1)
    
    if not type(maxambig) is float:
        try:
            maxambig = float(maxambig)
        except:
            sys.stderr.write(cyan(f"Error: --max-ambiguity must be between 0 and 1.\n"))
            sys.exit(-1)
    if not (maxambig <=1 and maxambig >=0):
        sys.stderr.write(cyan(f"Error: --max-ambiguity must be between 0 and 1.\n"))
        sys.exit(-1)

def fasta_ids_list(fasta):
    ids = []
    for record in SeqIO.parse(fasta,"fasta"):
        ids.append(record.id)
    ids = ids_qc(ids)
    return ids

def input_fasta_parsing(input_fasta,maxambig,minlen,config):
    
    misc.add_file_to_config(KEY_INPUT_SEQUENCES,input_fasta,config)
    misc.add_arg_to_config(KEY_MAX_AMBIGUITY,maxambig,config)
    misc.add_arg_to_config(KEY_MIN_LENGTH,minlen,config)
    
    if KEY_INPUT_SEQUENCES in config:
        if config[KEY_FROM_METADATA]:
            sys.stderr.write(cyan('Error: civet accepts either -fm/--from-metadata (which generates a query from the background data) or an input query (supplied via `-ids/--id-string`, `-i/--input-metadata` or `-f/--input-sequences`).\n'))
            sys.exit(-1)
        input_fasta_check(config[KEY_INPUT_SEQUENCES])
        fasta_qc_level(config[KEY_MAX_AMBIGUITY],config[KEY_MIN_LENGTH])

        if not KEY_INPUT_METADATA in config and not KEY_IDS in config:
            config[KEY_IDS] = fasta_ids_list(config[KEY_INPUT_SEQUENCES])

def parse_from_metadata(to_parse,background_csv):
    filters = {}
    column_names=''
    not_found = []

    if not type(to_parse)==list:
        to_parse = to_parse.split(" ")

    with open(background_csv,"r") as f:
        reader = misc.read_csv_or_tsv(background_csv,f)
        column_names = reader.fieldnames
        print(column_names)
        # get each of the factors for the query
        for factor in to_parse:
            # eg country=Ireland 
            column_name,to_search = factor.split("=")
            if column_name in column_names:
                filters[column_name] = to_search
            else:
                not_found.append(column_name)
                
    if len(not_found)==1:
        sys.stderr.write(cyan(f"Error: `-fm/--from-metadata` argument contains a column that is not found in the background metadata file supplied:") + f" {not_found[0]}\n")
        sys.exit(-1)
    elif len(not_found)>1:
        not_found_str = not_found.join("\n - ")
        sys.stderr.write(cyan(f"Error: `-fm/--from-metadata` argument contains columns not found in the background metadata file supplied:") + f"\n - {not_found_str}\n")
        sys.exit(-1)
    else:
        return filters

def check_date_format(date_string,row_number=None,column_name=None):
    date_format = '%Y-%m-%d'
    check_date= ""
    if date_string != "" and date_string != "NA":
        try:
            check_date = datetime.strptime(date_string, date_format).date()
        except:
            if row_number and column_name:
                sys.stderr.write(cyan(f"Error: Metadata field `{date_string}` [at column: {column_name}, row: {row_number}] contains unaccepted date format\nPlease use format {date_format}, i.e. `YYYY-MM-DD`\n"))
            else:
                sys.stderr.write(cyan(f"Error: Input '{date_string}' is the wrong date format.\nPlease use format {date_format}, i.e. `YYYY-MM-DD`\n"))

            sys.exit(-1)
            
    return check_date

def parse_date_range(background_csv,column_name,to_search,rows_to_search):
    date_range = to_search.split(":")
    start_date = datetime.strptime(date_range[0], "%Y-%m-%d").date()
    end_date = datetime.strptime(date_range[1], "%Y-%m-%d").date()

    if rows_to_search == []:
        with open(background_csv, newline="", encoding = "utf-8") as f:
            reader = misc.read_csv_or_tsv(background_csv,f)
            c =0
            for row in reader:
                c +=1
                row_date = row[column_name]
                
                check_date = check_date_format(row_date,c,column_name)

                if start_date <= check_date <= end_date:
                    rows_to_search.append((row,c))
    else:
        last_rows_to_search = rows_to_search
        new_rows_to_search = []
        for row,c in last_rows_to_search:
            row_date = row[column_name]

            check_date = check_date_format(row_date,c,column_name)

            if start_date <= check_date <= end_date:
                new_rows_to_search.append((row,c))

        rows_to_search = new_rows_to_search
    return rows_to_search

def parse_general_field(background_csv,column_name,to_search,rows_to_search):
    if rows_to_search == []:
        with open(background_csv, newline="", encoding="utf-8") as f:
            reader = misc.read_csv_or_tsv(background_csv,f)
            c =0
            for row in reader:
                c +=1
                row_info = row[column_name]
                
                if row_info.upper() == to_search:
                    rows_to_search.append((row,c))
    else:
        last_rows_to_search = rows_to_search

        new_rows_to_search = []
        for row,c in last_rows_to_search:
            row_info = row[column_name]
            
            if row_info.upper() == to_search:

                new_rows_to_search.append((row,c))

        rows_to_search = new_rows_to_search

    return rows_to_search

def filter_down_metadata(filters,background_csv):
    rows_to_search = []
    print(green("Query search:"))
    for column_name in filters:
        
        to_search = filters[column_name].upper()

        # assumes its a date range if it has a ':' and startswith 2020-, 2019- or 2021-
        if ':' in to_search:
            if to_search.startswith("2020-") or to_search.startswith("2019-") or to_search.startswith("2021-"):
                print(f"Date range: {to_search}")
                rows_to_search = parse_date_range(background_csv,column_name,to_search,rows_to_search)
        else:
            # parse by exact match 
            print(f"{column_name}: {to_search}")
            rows_to_search = parse_general_field(background_csv,column_name,to_search,rows_to_search)

    return rows_to_search

def from_metadata_parsing(config):
    if not config[KEY_FROM_METADATA] and not KEY_IDS in config:
        sys.stderr.write(cyan('Error: please specifiy a query or define a query using -fm/--from-metadata.\n'))
        sys.exit(-1)

    if config[KEY_FROM_METADATA]:
        filters = parse_from_metadata(config[KEY_FROM_METADATA],config["background_metadata"])
        # for column in [config["sequence_id_column"],config["input_id_column"],config["input_display_column"]]:
        #     column_names.append(column)

        rows_to_search = filter_down_metadata(filters,config["background_metadata"])
        
        query_ids = []
        # query_rows = []
        count = 0
        for row,c in rows_to_search:
            count +=1
            query_ids.append(row[config[KEY_BACKGROUND_ID_COLUMN]])
            # new_row = row
            # for column in [config["sequence_id_column"],config["input_id_column"],config["input_display_column"]]:
            #     new_row[colum] = row[config[KEY_BACKGROUND_ID_COLUMN]]
            # query_rows.append(new_row)

        if count > int(config[KEY_MAX_QUERIES]):
            sys.stderr.write(cyan(f'Error: -fm/--from-metadata found {count} matches, which exceeds the maximum query count.\n') + f"Either provide a more specific query with `-fm/--from-metadata` or overwrite the default maximum query limit with `-ql/--query-limit`.\n")
            sys.exit(-1)
        elif count == 0:
            sys.stderr.write(cyan(f"Error: No sequences meet the criteria defined with `--from-metadata`.\nPlease check your query is in the correct format (e.g. sample_date=YYYY-MM-DD).\nExiting\n"))
            sys.exit(-1)
        else:
            print(green(f"Number of sequences matching defined query:") + f" {count}")
        
        config[KEY_IDS] = query_ids

        # return query_ids,query_rows,column_names


