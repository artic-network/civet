import collections
import csv
from Bio import SeqIO
import hashlib

import os

def add_to_hash(seq_file,seq_map,hash_map,records):
    for record in SeqIO.parse(seq_file, "fasta"):
        records +=1
        seq = str(record.seq).encode()
        hash_object = hashlib.md5(seq)
        seq_map[hash_object.hexdigest()] = record.seq
        hash_map[hash_object.hexdigest()].append(record.id)
    return records

def merge(lists):
    newsets, sets = [set(lst) for lst in lists if lst], []
    while len(sets) != len(newsets):
        sets, newsets = newsets, []
        for aset in sets:
            for eachset in newsets:
                if not aset.isdisjoint(eachset):
                    eachset.update(aset)
                    break
            else:
                newsets.append(aset)
    return newsets

def get_merged_catchments(catchment_file,merged_catchment_file,config):
    
    query_dict = collections.defaultdict(list)
    
    # parse the catchment file
    with open(catchment_file, "r") as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            for col in ["closestsame","closestup","closestdown","closestside"]:
                seqs = row[col].split(";")
                for seq in seqs:
                    if seq:
                        query_dict[row["query"]].append(seq)
    
    # merge catchment sets
    lists = merge([query_dict[i] for i in query_dict])

    # organise the catchments and get out which query is in which catchment
    merged_catchments = {}
    catchment_key = {}
    catchment_dict = {}
    
    catchment_count = 0
    for i in lists:
        catchment_count +=1
        merged_catchments[f"catchment_{catchment_count}"] = i

        for query in query_dict:
            if query_dict[query][0] in i:
                catchment_key[query] = f"catchment_{catchment_count}"

    with open(merged_catchment_file,"w") as fw:
        fw.write(f'catchment,sequences\n')
        for catchment in merged_catchments:
            fw.write(f"{catchment},{';'.join(merged_catchments[catchment])}\n")

            for seq in merged_catchments[catchment]:
                catchment_dict[seq] = catchment

    return catchment_dict, catchment_key, catchment_count


def add_catchments_to_metadata(background_csv,query_metadata,query_metadata_with_catchments,catchment_dict,config):
    
    config["query_csv_header"].append("query_boolean")
    
    catchment_records = []

    with open(background_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            
            # exclude querys as they're already in the metadata
            if row[config["background_column"]] in config["ids"]:
                pass
            else:
                # the column used for sequence matching
                if row[config["fasta_column"]] in catchment_dict:
                    new_row = row

                    # query seqs excluded above
                    new_row["query_boolean"] = False

                    # add what catchment the sequence is in
                    new_row["catchment"] = catchment_dict[new_row[config["fasta_column"]]]

                    # check if it's got the headers from header in config, if not add them
                    for field in config["query_csv_header"]:
                        if field not in new_row:
                            new_row[field] = ""
                    catchment_records.append(new_row)

    with open(query_metadata_with_catchments,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=config["query_csv_header"],lineterminator='\n')
        writer.writeheader()

        with open(query_metadata, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                new_row = row
                new_row["query_boolean"] = True
                writer.writerow(new_row)

        for record in catchment_records:
            writer.writerow(record)

def write_catchment_fasta(catchment_dict,catchment_dir,config):
    seq_dict = collections.defaultdict(list)
    for record in SeqIO.parse(config["background_fasta"],"fasta"):
        if record.id in catchment_dict:
            seq_dict[catchment_dict[record.id]].append(record)

    for catchment in seq_dict:
        with open(os.path.join(catchment_dir,f"{catchment}.fasta"),"w") as fw:
            records = seq_dict[catchment]
            for record in SeqIO.parse(config["outgroup_fasta"],"fasta"):
                records.append(record)
            SeqIO.write(records,fw,"fasta")