#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
import sys
import os
import csv
import collections
from Bio import SeqIO
from time import time
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def account_for_all_ids(config):
    found_in_background_data = {}
    with open(config["background_metadata"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        for row in reader:
            if row[config["background_id_column"]] in config["ids"]:
                found_in_background_data[row[config["background_id_column"]]] = row
    
    print(green(f"Number of queries matched in background data:") + f" {len(found_in_background_data)}")
    
    found_in_input_fasta = {}
    if "input_sequences" in config:
        for record in SeqIO.parse(config["input_sequences"],"fasta"):
            if record.id in config["ids"]:
                found_in_input_fasta[record.id] = record
        print(green(f"Number of queries matched in supplied fasta:") + f" {len(found_in_input_fasta)}")
    
    not_found = []
    for i in config["ids"]:
        if i not in found_in_background_data and i not in found_in_input_fasta:
            not_found.append(i)

    if not_found:
        not_found_str = '\n - '.join(not_found)
        sys.stderr.write(cyan(f"Query records not matched:\n") + f"- {not_found_str}\n" + cyan("Please check input names match against records. Configure which columns to match using `-icol/--input-id-column` and `-bicol/--background-id-column`.\n"))
        sys.exit(-1)
    
    return found_in_background_data,header,found_in_input_fasta

def merge_metadata_records(found_in_background_data,header,config):
    record_info = collections.defaultdict(dict)

    for record in found_in_background_data:
        record_info[record] = found_in_background_data[record]
        if not config["input_id_column"] in record_info[record]:
            record_info[record][config["input_id_column"]] = record_info[record][config["background_id_column"]]

    if "input_metadata" in config:
        with open(config["input_metadata"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row[config["input_id_column"]]
                for col in [config["background_id_column"],config["sequence_id_column"]]:
                    if not col in row:
                        record_info[name][col] = name

                for key in row:
                    if key not in header:
                        header.append(key)
                    
                    record_info[name][key] = row[key]
    else:
        for record in config["ids"]:
            for col in [config["input_id_column"],config["background_id_column"]]: 
                if col not in header:
                    header.append(col)
                record_info[record][col] = record
    
    for key in header:
        for record in record_info:

            if key not in record_info[record]:
                record_info[record][key]=""
    
    return record_info,header

def add_sequence_status_to_metadata(record_info,header,found_in_background_data,passed_qc,failed_qc):
    header.append("source")
    header.append("qc_status")
    header.append("seq_N_content")
    header.append("seq_length")

    # creating transiently the same data structure as the failed_qc seqs
    passed_info = {}
    for record in passed_qc:
        passed_info[record[0].id] = (record[0].id, record[1],record[2])

    for record in record_info:
        if record in found_in_background_data:
            record_info[record]["source"] = "background_data"
            record_info[record]["qc_status"] = ""
            record_info[record]["seq_N_content"] = ""
            record_info[record]["seq_length"] = ""
        else:
            record_info[record]["source"] = "input_fasta"

            if record in failed_qc:
                record_info[record]["qc_status"] = f"Fail {failed_qc[record][0]}"
                record_info[record]["seq_N_content"] = failed_qc[record][1]
                record_info[record]["seq_length"] = failed_qc[record][2]
            else:
                record_info[record]["qc_status"] = "Pass"
                record_info[record]["seq_N_content"] = passed_info[record][1]
                record_info[record]["seq_length"] = passed_info[record][2]

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
            
            passed_ids = [i[0].id for i in passed_qc]

            num_N = str(record.seq).upper().count("N")
            prop_N = round((num_N)/len(record.seq), 2)
            if len(record) < minlen:
                # fail[id] = (why,n content,len)
                failed_qc[record.id] = (f"Sequence too short: Sequence length {len(record)}",prop_N,len(record))
            else:
                if prop_N > maxambig: 
                    failed_qc[record.id] = (f"N content: {prop_N}",prop_N,len(record))
                else:
                    if record.id in passed_ids:
                        failed_qc[record.id] = ("Sequence duplicated in input fasta file.",prop_N,len(record))
                    else:
                        passed_qc.append((record,prop_N,len(record)))
        if len(passed_qc) == 1:
            print(green(f"{len(passed_qc)} supplied query sequence has passed QC."))
        else:
            print(green(f"{len(passed_qc)} supplied query sequences have passed QC."))

        if len(failed_qc) == 1:
            print(cyan(f"{len(failed_qc)} supplied query sequence has failed QC."))
        else:
            print(cyan(f"{len(failed_qc)} supplied query sequences have failed QC."))

        if len(passed_qc) == 0:
            sys.stderr.write(cyan(f"No supplied sequences pass the QC steps.\n") + f"""Please ensure sequences meet the following:
    \t- Minimum sequence length (>{minlen} bases)
    \t- Maximum ambiguities (<{maxambig} proportion N)
    \t- Whether the records match a query supplied in an ID string or input csv
    \t- Whether the records are duplicated in the file\n""" + cyan("You can change the default QC settings with `-n/--max-ambiguity` and `-l/--min-length`."))
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

    return query_metadata, passed_qc, found_in_background_data

def write_master_metadata(query_metadata, config):
    """
    Requires directory_setup to have been run
    """
    
    with open(os.path.join(config["data_outdir"],"query_metadata.master.csv"),"w") as fw:
        
        writer = csv.DictWriter(fw, fieldnames=config["query_csv_header"],lineterminator='\n')
        writer.writeheader()

        for row in query_metadata:
            writer.writerow(query_metadata[row])

    config["query_metadata"] = os.path.join(config["data_outdir"],"query_metadata.master.csv")

def write_passed_qc_fasta(passed_qc, config):
    if passed_qc:
        passed_records = [i[0] for i in passed_qc]
        config["query_fasta"] = os.path.join(config["tempdir"],"query.passed_qc.fasta")
        with open(config["query_fasta"],"w") as fw:
            SeqIO.write(passed_records, fw, "fasta")
    else:
        config["query_fasta"] = False

def write_matched_fasta(found_in_background_data, config):

    num_found_in_background_data = len(found_in_background_data)
    sequence_names_to_match = [found_in_background_data[i][config["sequence_id_column"]] for i in found_in_background_data]
    sequence_names_to_match = set(sequence_names_to_match)
    config["matched_fasta"] = False
    if num_found_in_background_data != 0:
        matched_records = []

        start = time()
        for record in SeqIO.parse(config["background_sequences"],"fasta"):
            if len(matched_records)!= num_found_in_background_data:
                if record.id in sequence_names_to_match:
                    matched_records.append(record)
            else:
                break

        end = time()
        # print("Biopython")
        # print("Start:",start,"End:",end, "Diff:",end-start)
        # start = time()
        # matched_records = []
        # with open(config["background_sequences"],"r") as f:
        #     record_id = ""
        #     record_seq = ""
        #     new_id = ""
        #     for l in f:
        #         l=l.rstrip("\n")
        #         if len(matched_records)!= num_found_in_background_data:
        #             if l[0]=='>':
        #                 name=l.lstrip(">")
        #                 if name in sequence_names_to_match:
        #                     record_id = new_id
        #                     new_id =f"{name}"
        #                     record_seq = ""
        #                     continue
        #                 else:
        #                     continue

        #             if not l[0]=='>':
        #                 record_seq+=l
                    
        #             record = SeqRecord(Seq(record_seq), id=record_id)
        #             matched_records.append(record)
        #         else:
        #             break
                
        # end = time()
        # print("Base")
        # print("Start:",start,"End:",end, "Diff:",end-start)

        if not matched_records:
            sys.stderr.write(cyan(f"""Error: No sequence records matched.\nPlease check the `-sicol/--sequence-id-column` is matching the sequence ids.
Currently searching the `{config['sequence_id_column']}` column.\n"""))
            sys.exit(-1)
        else:
            if len(matched_records) == num_found_in_background_data:
                config["matched_fasta"] = os.path.join(config["tempdir"],"query.matched_background.fasta")
                with open(os.path.join(config["matched_fasta"]),"w") as fw:
                    SeqIO.write(matched_records, fw, "fasta")
            else:
                not_found = []
                matched_ids = [record.id for record in matched_records]
                for record in found_in_background_data:
                    if record not in matched_ids:
                        not_found.append(record)
                not_found_str = "\n- ".join(not_found)
                sys.stderr.write(cyan(f"Some records in background metadata file not found in background fasta file.\n") + f"- {not_found_str}\n" + cyan("Please check sequence ids match against the data in `-bicol/--background-id-column`.\n"))
                sys.exit(-1)

def write_parsed_query_files(query_metadata,passed_qc,found_in_background_data, config):
    write_master_metadata(query_metadata, config)

    write_passed_qc_fasta(passed_qc, config)

    write_matched_fasta(found_in_background_data, config)