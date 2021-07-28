import collections
import csv
from Bio import SeqIO
import hashlib
from numpy.random import choice
import sys
from civet.utils.log_colours import green,cyan,red
import os
csv.field_size_limit(sys.maxsize)

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
            query_dict[row["query"]].append(row["query"])
            for col in ["closestsame","closestup","closestdown","closestside"]:
                seqs = row[col].split(";")
                for seq in seqs:
                    if seq:
                        if seq not in config["ids"]:
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

        if len(i) == 1:
            sys.stderr.write(cyan(f'catchment_{catchment_count} has only one item, please increase the SNP distance to get larger catchments.\n'))
            sys.exit(0)

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
            if row[config["background_id_column"]] in config["ids"]:
                pass
            else:
                # the column used for sequence matching
                if row[config["sequence_id_column"]] in catchment_dict:
                    new_row = row

                    # query seqs excluded above
                    new_row["query_boolean"] = "False"

                    # add what catchment the sequence is in
                    new_row["catchment"] = catchment_dict[new_row[config["sequence_id_column"]]]
                    new_row[config["input_display_column"]] = new_row[config["sequence_id_column"]]

                    # check if it's got the headers from header in config, if not add them
                    for field in config["query_csv_header"]:
                        if field not in new_row:
                            new_row[field] = ""
                    for mutation in config["mutations"]:
                        if mutation not in new_row:
                            new_row[mutation] = "False"
                    catchment_records.append(new_row)

    with open(query_metadata_with_catchments,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=config["query_csv_header"],lineterminator='\n')
        writer.writeheader()

        with open(query_metadata, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                new_row = row
                new_row["query_boolean"] = "True"
                writer.writerow(new_row)

        for record in catchment_records:
            writer.writerow(record)

def smooth_weights(column,metadata):
    counter = collections.Counter()
    for row in metadata:
        counter[row[column]]+=1
    prop = 1/len(counter)
    weights = {}
    for item in counter:
        weights[item]=prop/counter[item]
    weight_list = []
    for row in metadata:
        weight_list.append(weights[row[column]])
    return weight_list

def enrich_weights(column,field,factor,metadata):
    cum_weight = 0
    row_dict = {}
    for row in metadata:
        if row[column] == field:
            cum_weight +=1*factor
            row_dict[row[column]] = 1*factor
        else:
            cum_weight +=1
            row_dict[row[column]] = 1
        
    weights = {}
    for item in row_dict:
        weights[item]=row_dict[item]/cum_weight
    
    weight_list = []
    for row in metadata:
        weight_list.append(weights[row[column]])
    return weight_list

def downsample_catchment(catchment_metadata,size,strategy,column,background_column,field=None,factor=None):
    names = [row[background_column] for row in catchment_metadata]
    if strategy == "random":
        downsample = choice(names,size=size,replace=False)
    else:
        if strategy == "normalise":
            weights = smooth_weights(column,catchment_metadata)
        elif strategy == "enrich":
            weights = enrich_weights(column,field,factor,catchment_metadata)
        
        downsample = choice(names,size=size,replace=False,p=weights)

    downsample_metadata = []
    for row in catchment_metadata:
        new_row = row
        if row[background_column] in downsample:
            new_row["in_tree"] = "True"
        else:
            new_row["in_tree"] = "False"
        downsample_metadata.append(new_row)

    return downsample_metadata

def write_catchment_fasta(catchment_metadata,fasta,catchment_dir,config):
    catchment_dict = {}
    
    with open(catchment_metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["catchment"] in config["figure_catchments"]:
                if row["query_boolean"] == "True":
                    catchment_dict[row["hash"]] = row["catchment"]
                else:
                    catchment_dict[row[config["sequence_id_column"]]] = row["catchment"]

    seq_dict = collections.defaultdict(list)
    seq_count = 0
    for record in SeqIO.parse(fasta,"fasta"):
        if record.id in catchment_dict:
            seq_count+=1
            seq_dict[catchment_dict[record.id]].append(record)

    for record in SeqIO.parse(config["background_sequences"],"fasta"):
        if seq_count != len(catchment_dict):
            if record.id in catchment_dict:
                seq_dict[catchment_dict[record.id]].append(record)
                seq_count +=1
        else:
            break

    if seq_count == 0:
        sys.stderr.write(cyan(f"Error: No sequence records matched.\nPlease check the `-sicol/--sequence-index-column` is matching the sequence ids.\n"))
        sys.exit(-1)
    
    for catchment in seq_dict:

        with open(os.path.join(catchment_dir,f"{catchment}.fasta"),"w") as fw:
            records = seq_dict[catchment]
            for record in SeqIO.parse(config["outgroup_fasta"],"fasta"):
                records.append(record)
            SeqIO.write(records,fw,"fasta")

def which_catchments_too_large(in_csv,config):
    catchment_counter= collections.Counter()

    figure_catchments = set()
    with open(in_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["query_boolean"] == "True":
                catchment_counter[row["catchment"]] +=1


    with open(in_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            catchment  = row["catchment"]
            if catchment_counter[catchment] <= int(config["max_tree_size"]):
                figure_catchments.add(catchment)

    config["figure_catchments"] = list(figure_catchments)


def downsample_if_building_trees(in_csv, out_csv, config):
    with open(out_csv,"w") as fw:
        with open(in_csv,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            # in tree to header
            header.append("in_tree")

            writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
            writer.writeheader()
            
            if not '3' in config["report_content"]:
                # if no tree, don't need to downsample- just write all true
                for row in reader:
                    new_row= row
                    new_row["in_tree"]="False"
                    writer.writerow(new_row)
            else:
                # if going to build tree, see if need to downsample
                catchment_dict = collections.defaultdict(list)
                
                for row in reader:
                    if row["catchment"] in config["figure_catchments"]:
                        if row["query_boolean"] == "True":
                            new_row = row
                            new_row["in_tree"]="True"
                            # write out the query metadata to file
                            writer.writerow(new_row)
                        else:
                            # categorise the metadata by catchment
                            catchment_dict[row["catchment"]].append(row)
                    else:
                        new_row = row
                        new_row["in_tree"]="False"
                        writer.writerow(new_row)
                
                for catchment in catchment_dict:

                    if catchment not in config["figure_catchments"]:
                        # dont build tree
                        downsample_metadata = []
                        for row in catchment_dict[catchment]:
                            new_row = row
                            new_row["in_tree"] = "False"
                            downsample_metadata.append(new_row)
                            #write out the new info for non queries
                    else:                        
                        downsample = False
                        # figure out how many seqs to downsample to per catchment
                        target = config["catchment_background_size"]

                        # if len(catchment_dict[catchment]) > config["catchment_background_size"]:
                        #     print(cyan(f'Warning: number of queries in {catchment} (n={query_counter[catchment]}) exceeds the maximum catchment size of {config["catchment_background_size"]}.\nPlease be aware tree building may take longer than expected.'))
                            
                        if len(catchment_dict[catchment]) > target:
                            downsample = True

                        if downsample == True:
                            # need to downsample
                            downsample_metadata = downsample_catchment(catchment_dict[catchment],target,config["mode"],config["downsample_column"],config["sequence_id_column"],config["downsample_field"],config["factor"])
                        else:
                            # dont need to downsample
                            downsample_metadata = []
                            for row in catchment_dict[catchment]:
                                new_row = row
                                new_row["in_tree"] = "True"
                                downsample_metadata.append(new_row)
                                #write out the new info for non queries

                    for new_row in downsample_metadata:
                        writer.writerow(new_row)


def which_trees_to_run(in_csv):
    catchment_counter= collections.Counter()
    with open(in_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["in_tree"] == "True":
                catchment_counter[row["catchment"]] +=1
    to_run = []
    for catchment in catchment_counter:
        if catchment_counter[catchment] != 0:
            to_run.append(catchment)
    return ','.join(to_run)