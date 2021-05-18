#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
import sys
import os
import csv
import collections
from Bio import SeqIO

def account_for_all_ids(config):
    found_in_background_data = {}
    with open(config["background_csv"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        for row in reader:
            if row[config["data_column"]] in config["ids"]:
                found_in_background_data[row[config["data_column"]]] = row
    
    print(green(f"Number of queries matched in background data:") + f" {len(found_in_background_data)}")
    
    found_in_input_fasta = {}
    if "fasta" in config:
        for record in SeqIO.parse(config["fasta"],"fasta"):
            if record.id in config["ids"]:
                found_in_input_fasta[record.id] = record
    print(green(f"Number of queries matched in supplied fasta:") + f" {len(found_in_input_fasta)}")
    
    not_found = []
    for i in config["ids"]:
        if i not in found_in_background_data and i not in found_in_input_fasta:
            not_found.append(i)

    if not_found:
        not_found_str = '\n - '.join(not_found)
        sys.stderr.write(cyan(f"Query records not matched:\n") + f"- {not_found_str}\n" + cyan("Please check input names match against records. Configure which columns to match using `-icol/--input-column` and `-dcol/--data-column`.\n"))
        sys.exit(-1)

    return found_in_background_data,header,found_in_input_fasta

def merge_metadata_records(found_in_background_data,header,config):
    record_info = collections.defaultdict(dict)

    for record in found_in_background_data:
        record_info[record] = found_in_background_data[record]

    if "input_csv" in config:
        with open(config["input_csv"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row[config["input_column"]]
                for key in row:
                    if key not in header:
                        header.append(key)
                    record_info[name][key] = row[key]
    else:
        for record in config["ids"]:
            record_info[record][config["input_column"]] = record
    
    for key in header:
        for record in record_info:
            if key not in record_info[record]:
                record_info[record][key]=""
    
    return record_info,header

def add_sequence_status_to_metadata(record_info,header,found_in_background_data,passed_qc,failed_qc):
    header.append("source")
    header.append("qc_status")
    for record in record_info:
        if record in found_in_background_data:
            record_info[record]["source"] = "background_data"
            record_info[record]["qc_status"] = ""
        else:
            record_info[record]["source"] = "input_fasta"

            if record in failed_qc:
                record_info[record]["qc_status"] = f"Fail {failed_qc[record]}"
            else:
                record_info[record]["qc_status"] = "Pass"

    return record_info, header

def input_fasta_qc(found_in_input_fasta,config):
    """
    Checking input fasta file for:
    - Minimum sequence length
    - Maximum ambiguities
    - Whether the records match a query id
    - Whether the records are duplicated
    """
    failed_qc = {}
    passed_qc = []

    if found_in_input_fasta:
        minlen = config["min_length"]
        maxambig = config["max_ambiguity"]
        for record_name in found_in_input_fasta:
            record = found_in_input_fasta[record_name]
            if len(record) < minlen:
                failed_qc[record.id] = f"Sequence too short: Sequence length {len(record)}"
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > maxambig: 
                    failed_qc[record.id] = f"N content: {prop_N}"
                else:
                    if record.id in passed_qc:
                        failed_qc[record.id] = "Sequence duplicated in input fasta file."
                    else:
                        passed_qc.append(record)
        
        print(green(f"{len(passed_qc)} supplied query sequences have passed QC."))
        print(cyan(f"{len(failed_qc)} supplied query sequences have failed QC."))

        if len(passed_qc) == 0:
            sys.stderr.write(cyan(f"No supplied sequences pass the QC steps.\n") + f"""Please ensure sequences meet the following:
    \t- Minimum sequence length (>{minlen} bases)
    \t- Maximum ambiguities (<{maxambig} proportion N)
    \t- Whether the records match a query supplied in an ID string or input csv
    \t- Whether the records are duplicated in the file\n""")
            sys.exit(-1)

    return passed_qc, failed_qc

def check_if_any_valid_queries(found_in_background_data,passed_qc):

    if len(found_in_background_data) + len(passed_qc) == 0:
        sys.stderr.write(cyan(f"No query records identified or passed QC\nExiting."))
        sys.exit(-1)

def query_check_against_background_merge_input(config):
    # how many ids input
    # match against background
    # match against background fasta
    
    found_in_background_data,header,found_in_input_fasta = account_for_all_ids(config)

    # supplied fasta qc
    passed_qc,failed_qc = input_fasta_qc(found_in_input_fasta,config)

    # if nothing in background matched and nothing in supplied fasta matched, exit
    check_if_any_valid_queries(found_in_background_data,passed_qc)

    # merge metadata from background and input, overwrite with supplied if overlap
    # record_dict[input_col_data] = {'col1':'data','col2':'data'...}
    query_metadata,header = merge_metadata_records(found_in_background_data,header,config)

    query_metadata, header = add_sequence_status_to_metadata(query_metadata,header,found_in_background_data,passed_qc,failed_qc)
    config["query_csv_header"] = header

    return query_metadata, passed_qc
