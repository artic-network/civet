#!/usr/bin/env python3

import argparse
import collections
from Bio import SeqIO
import os
import csv
cwd = os.getcwd()



def parse_args():
    parser = argparse.ArgumentParser(description='Check COG DB for query sequences.')

    parser.add_argument("--query", action="store", type=str, dest="query")
    parser.add_argument("--cog-seqs", action="store", type=str, dest="cog_seqs")
    parser.add_argument("--cog-metadata", action="store", type=str, dest="cog_metadata")
    parser.add_argument("--field", action="store", type=str, dest="field")
    parser.add_argument("--in-metadata", action="store", type=str, dest="in_metadata")
    parser.add_argument("--in-seqs", action="store", type=str, dest="in_seqs")
    parser.add_argument("--not-in-cog", action="store", type=str, dest="not_in_cog")
    parser.add_argument("--all-cog",action="store_true",dest="all_cog")
    return parser.parse_args()

def check_cog_db():
    args = parse_args()

    found = []
    
    in_cog_metadata = []
    in_cog_names = {}

    column_to_match = args.field
    
    query_names = []
    with open(args.query,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            query_names.append(row["name"])

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
                    if args.all_cog:
                        row["source"]="COG database"
                    else:
                        row["source"]="phylogeny"
                    in_cog_metadata.append(row)
                    in_cog_names[row[column_to_match]] = row["sequence_name"]
    if args.all_cog:
        print(f"Number of query records found in CLIMB: {len(in_cog_metadata)}")
    else:
        print(f"Number of query records found in tree: {len(in_cog_metadata)}")
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
    print(f"Number of associated sequences found: {len(found)}")
    with open(args.not_in_cog, "w") as fw:
        if args.all_cog:
            print("\nThe following sequences were not found in the cog database:")
        else:
            print("\nThe following sequences were not found in the phylogeny:")
        
        fw.write("name\n")
        for query in query_names:
            if query not in found:
                fw.write(query + '\n')
                print(f"{query}")
        print("If you wish to access sequences in the cog database\nwith your query, ensure you have the correct sequence id.")


if __name__ == '__main__':

    check_cog_db()