#!/usr/bin/env python3
import log_colours as colour
import sys
import os

def make_csv_string_from_ids(ids, config):
    ids = ids.split(",")
    with open(query,"w") as fw:
        in_col = "name"
        config["input_column"] = in_col
        fw.write(f"{in_col}\n")
        c = 0
        for i in id_list:
            c +=1
            fw.write(i+'\n')
        print(green(f"Number of IDs:") + f" {c}")
    return query

def input_csv_qc(cwd,inputfile):
    input_file = os.path.join(cwd,inputfile)

    ending = input_file.split(".")[-1]

    if ending == "csv":
        print(colour.green(f"Input file:") + f" {input_file}")
        query = input_file

    elif ending == "xls":
        sys.stderr.write(colour.cyan(f"Error: it looks like you've provided an excel file as input. Please provide a csv file.\n"))
        sys.exit(-1)
    else:
        sys.stderr.write(colour.cyan(f"Error: -i,--input accepts a csv file, please provide this or a comma-separated string of IDs.\n"))
        sys.exit(-1)

