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
