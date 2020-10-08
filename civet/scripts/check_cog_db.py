#!/usr/bin/env python3

import argparse
import collections
from Bio import SeqIO
import os
import csv
cwd = os.getcwd()

from reportfunk.funks import io_functions as qcfunk

def parse_args():
    parser = argparse.ArgumentParser(description='Check COG DB for query sequences.')

    parser.add_argument("--query", action="store", type=str, dest="query")
    parser.add_argument("--cog-seqs", action="store", type=str, dest="cog_seqs")
    parser.add_argument("--cog-metadata", action="store", type=str, dest="cog_metadata")
    parser.add_argument("--field", action="store", type=str, dest="field")
    parser.add_argument("--in-metadata", action="store", type=str, dest="in_metadata")
    parser.add_argument("--in-seqs", action="store", type=str, dest="in_seqs")
    parser.add_argument("--not-in-cog", action="store", type=str, dest="not_in_cog")
    parser.add_argument("--input-column",action="store",dest="input_column")
    return parser.parse_args()

def check_cog_db():
    args = parse_args()

    found = []
    
    in_cog_metadata = []
    in_cog_names = {}

    column_to_match = args.field
    input_column = args.input_column
    query_names = []
    with open(args.query,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            query_names.append(row[input_column])

    with open(args.cog_metadata,newline="") as f:
        reader = csv.DictReader(f)
        header_names = reader.fieldnames
        for row in reader:
            for seq in query_names:

                cog_id = row[column_to_match]
                if seq == cog_id:
                    
                    row["query_id"]=seq
                    row["cog_id"] = row[column_to_match]
                    row["query"]=row["sequence_name"]
                    row["closest"]=row["sequence_name"]
                    row["source"]="phylogeny"
                    in_cog_metadata.append(row)
                    in_cog_names[row[column_to_match]] = row["sequence_name"]
    
    print(qcfunk.green(f"Number of query records found in tree:")+f" {len(in_cog_metadata)}")
    with open(args.in_metadata, "w") as fw:
        header_names.append("query_id")
        header_names.append("cog_id")
        header_names.append("query")
        header_names.append("closest")
        header_names.append("source")
        writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
        writer.writeheader()
        writer.writerows(in_cog_metadata)

    with open(args.in_seqs, "w") as fw:
        for record in SeqIO.parse(args.cog_seqs, "fasta"):
            for name in in_cog_names:
                sequence_name = in_cog_names[name]
                if sequence_name==record.id:
                    found.append(name)
                    fw.write(f">{name} sequence_name={record.id} status=in_cog\n{record.seq}\n")

    with open(args.not_in_cog, "w") as fw:
        
        c = 0
        not_found_str = ""
        fw.write(f"{input_column}\n")
        for query in query_names:
            if query not in found:
                fw.write(query + '\n')
                not_found_str += (f"\t- {query}\n")
                c+=1
        if c != 0:
            print(qcfunk.cyan("\nNot found in phylogeny:"))
            print(not_found_str)



if __name__ == '__main__':

    check_cog_db()