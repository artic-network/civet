#!/usr/bin/env python3
import csv
import argparse
from Bio import SeqIO
import collections
"""
--paf {input.paf:q} 
--metadata {input.metadata:q} 
--csv-out {output.csv:q} 
"""
def parse_args():
    parser = argparse.ArgumentParser(description='Parse barcode info and minimap paf file, create report.')

    parser.add_argument("--paf", action="store", type=str, dest="paf")
    parser.add_argument("--metadata", action="store", type=str, dest="metadata")
    parser.add_argument("--search-field", action="store",type=str, dest="search_field")
    parser.add_argument("--csv-out", action="store", type=str, dest="outfile")
    parser.add_argument("--seqs", action="store", type=str, dest="seqs")
    parser.add_argument("--seqs-out", action="store", type=str, dest="seqs_out")
    return parser.parse_args()


def parse_line(line):
    values = {}
    
    tokens = line.rstrip('\n').split('\t')
    values["name"], values["read_len"] = tokens[:2]
    values["query_start"] = tokens[2]
    values["query_end"] = tokens[3]
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"] = tokens[5:11]

    return values


def get_closest_cog_sequences(paf):

    closest_cog_sequences = []
    closest_to_query = collections.defaultdict(list)
    with open(paf,"r") as f:
        last_mapping = None
        for line in f:

            mapping = parse_line(line)

            if last_mapping:
                if mapping["name"] == last_mapping["name"]:
                    last_mapping['ref_hit'] = '?'
                else:
                    closest_cog_sequences.append(last_mapping["ref_hit"])
                    closest_to_query[last_mapping["ref_hit"]].append(last_mapping["name"])
                    last_mapping = mapping
            else:
                last_mapping = mapping

        closest_cog_sequences.append(last_mapping["ref_hit"])
        closest_to_query[last_mapping["ref_hit"]].append(last_mapping["name"])

    return closest_cog_sequences, closest_to_query


def parse_paf_and_get_metadata():
    args = parse_args()

    closest_cog, closest_to_query = get_closest_cog_sequences(args.paf)
    column_to_match = args.search_field
    with open(args.metadata, newline="") as f:
        rows_to_write = []
        reader = csv.DictReader(f)
        header_names = reader.fieldnames
        with open(args.outfile, "w") as fw:
            header_names.append("query_id")
            header_names.append("cog_id")
            header_names.append("query")
            header_names.append("closest")
            writer = csv.DictWriter(fw, fieldnames=header_names)
            writer.writeheader()
        
            for row in reader:
                if row["sequence_name"] in closest_cog:
                    for query in closest_to_query[row["sequence_name"]]:
                        new_row = row
                        new_row["query_id"]=query
                        new_row["cog_id"]= row[column_to_match]
                        new_row["query"]=query
                        new_row["closest"]=row["sequence_name"]

                        writer.writerow(new_row)


    with open(args.seqs_out, "w") as fw:
        for record in SeqIO.parse(args.seqs,"fasta"):
            if record.id in closest_cog:
                closest_queries = ",".join(closest_to_query[record.id])
                fw.write(f">{record.id} query={closest_queries}\n{record.seq}\n")

if __name__ == '__main__':

    parse_paf_and_get_metadata()
    