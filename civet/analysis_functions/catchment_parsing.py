import collections
import csv
from Bio import SeqIO
import hashlib


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

def get_merged_catchments(catchment_file,key_file,merged_catchment_file,config):
    
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
    
    catchment_count = 0
    for i in lists:
        catchment_count +=1
        merged_catchments[f"catchment_{catchment_count}"] = i

        for query in query_dict:
            if query_dict[query][0] in i:
                catchment_key[query] = f"catchment_{catchment_count}"
    
    with open(key_file,"w") as fw:
        fw.write(f'{config["input_column"]},catchment\n')
        for query in catchment_key:
            fw.write(f"{query},{catchment_key[query]}\n")

    with open(merged_catchment_file,"w") as fw:
        fw.write(f'catchment,sequences\n')
        for catchment in merged_catchments:
            fw.write(f"{catchment},{';'.join(merged_catchments[catchment])}\n")

    return merged_catchments

