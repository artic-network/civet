#!/usr/bin/env python3
from collections import defaultdict
import pandas as pd
import csv
from tabulate import tabulate
import baltic as bt
import os
import datetime as dt


class taxon():

    def __init__(self, name, global_lin, uk_lin, phylotype):

        self.name = name

        self.sample_date = "NA"
        
        if global_lin == "":
            self.global_lin = "NA"
        else:
            self.global_lin = global_lin
        
        if uk_lin == "":
            self.uk_lin = "NA"
        else:
            self.uk_lin = uk_lin
       
        if phylotype == "":
            self.phylotype = "NA"
        else:
            self.phylotype = phylotype
       
        self.in_cog = False

        self.attribute_dict = {}


def parse_reduced_metadata(metadata_file):
    
    query_dict = {}
    query_id_dict = {}
    present_lins = set()

    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}

    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            glob_lin = sequence['lineage']
            uk_lineage = sequence['uk_lineage']
            
            adm1_prep = sequence["adm1"].split("-")[1]
            adm1 = contract_dict[adm1_prep]

            query_id = sequence['query_id']
            query_name = sequence['query']
            closest_name = sequence["closest"]

            phylotype = sequence["phylotype"]
            sample_date = sequence["sample_date"]

            new_taxon = taxon(query_name, glob_lin, uk_lineage, phylotype)

            new_taxon.query_id = query_id

            present_lins.add(uk_lineage)
            

            if query_name == closest_name: #if it's in COG, get it's sample date
                new_taxon.in_cog = True
                #new_taxon.attribute_dict["sample_date"] = sample_date
                new_taxon.sample_date = sample_date
                new_taxon.closest = "NA"
            else:
                new_taxon.closest = closest_name

            new_taxon.attribute_dict["country"] = "UK"
            new_taxon.attribute_dict["adm1"] = adm1

            query_dict[query_name] = new_taxon
            query_id_dict[query_id] = new_taxon
            
    return query_dict, query_id_dict, present_lins

def parse_input_csv(input_csv, query_id_dict, desired_fields):

    new_query_dict = {}
    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIRE": "Northern_Ireland"}

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            name = sequence["name"]

            if name in query_id_dict.keys():
                taxon = query_id_dict[name]

                if "sample_date" in col_names:
                    #taxon.attribute_dict["sample_date"] = sequence["sample_date"] #if it's not in COG but date is provided
                    if sequence["sample_date"] != "":
                        taxon.sample_date = sequence["sample_date"]
                
                #if it's in COG, it will already have been assigned a sample date.

                for col in col_names:
                    if desired_fields != []:
                        if col != "name" and col in desired_fields:
                            if sequence[col] == "":
                                taxon.attribute_dict[col] = "NA"
                            else:
                                taxon.attribute_dict[col] = sequence[col]
            
                        elif col == "adm1":
                            if "UK" in sequence[col]:
                                adm1_prep = sequence[col].split("-")[1]
                                adm1 = contract_dict[adm1_prep]
                            else:
                                adm1 = sequence[col]

                            taxon.attribute_dict["adm1"] = adm1

                new_query_dict[taxon.name] = taxon
            
            # else:
            #     print(name + " is in the input file but not the processed file. This suggests that it is not in COG and a sequence has also not been provided.")
                
    return new_query_dict 


def parse_tree_tips(tree_dir):

    tips = []

    for fn in os.listdir(tree_dir):
        if fn.endswith("tree"):
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
            for k in tree.Objects:
                if k.branchType == 'leaf' and "inserted" not in k.name:
                    tips.append(k.name)

        elif fn.endswith(".txt"):
            with open(tree_dir + "/" + fn) as f:
                for l in f:
                    tip_string = l.strip("\n").split("\t")[1]
                    tip_list = tip_string.split(",")
                    tips.extend(tip_list)

    return tips

def parse_full_metadata(query_dict, full_metadata, present_lins, present_in_tree):

    full_tax_dict = query_dict.copy()

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            uk_lin = sequence["uk_lineage"]
            seq_name = sequence["sequence_name"]

            date = sequence["sample_date"]
            adm2 = sequence["adm2"]
            country = sequence["country"]

            glob_lin = sequence["lineage"]
            phylotype = sequence["phylotype"]

            if (uk_lin in present_lins or seq_name in present_in_tree) and seq_name not in query_dict.keys():
                new_taxon = taxon(seq_name, glob_lin, uk_lin, phylotype)
                if date == "":
                    date = "NA"
                
                
                # new_taxon.attribute_dict["sample_date"] = date
                new_taxon.sample_date = date
                
                new_taxon.attribute_dict["adm2"] = adm2
                new_taxon.attribute_dict["country"] = country

                full_tax_dict[seq_name] = new_taxon

    return full_tax_dict

def make_initial_table(query_dict):

    df_dict = defaultdict(list)

    for query in query_dict.values():
        
        df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        df_dict["Part of COG"].append(query.in_cog)
        
        if query.in_cog: 
            df_dict["Sequence name in COG"].append(query.name)
        else:
            df_dict["Sequence name in COG"].append("NA")

        df_dict["Closest sequence in COG"].append(query.closest)
        
        df_dict["UK lineage"].append(query.uk_lin)
        df_dict["Global lineage"].append(query.global_lin)
        df_dict["Phylotype"].append(query.phylotype)

        for key, value in query.attribute_dict.items():
            df_dict[key].append(value)

    df = pd.DataFrame(df_dict)

    df.set_index("Query ID", inplace=True)

    return df

def investigate_QC_fails(QC_file):

    fail_dict = {}

    with open(QC_file) as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            name = sequence["name"]
            reason = sequence["reason_for_failure"]

            if "seq_len" in reason:
                length = reason.split(":")[1]
                final_reason = "Sequence too short: only " + length + " bases."
            elif "N_content" in reason:
                n_content = reason.split(":")[1]
                final_reason = "Sequence has too many Ns: " + str(float(n_content)*100) + "\% of bases"

            fail_dict[name] = final_reason


    return fail_dict

def print_missing_seqs(missing_seqs_file):
    
    failed_names = []

    with open(missing_seqs_file) as f:
        for l in f:
            name = l.strip("\n").split(",")[0]
            
            failed_names.append(name)

    return failed_names


